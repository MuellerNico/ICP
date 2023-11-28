using Distributed


function init_grid!(grid, origin, Δx::Float64, A::Float64, sigma::Float64)
    for j in 1:size(grid)[2]
        # -1: The top left corner of the grid is at the origin.
        # However, the index of the first row is 1.
        x = origin[1] + Δx * (j-1)
        for i in 1:size(grid)[1]
            y = origin[2] + Δx * (i-1)
            
            r_squared = x^2 + y^2
            grid[i, j, 1] = A * exp(-r_squared / sigma^2)
        end
    end
end


function calculate_RHS_ij(grid, Δx::Float64, i::Int64, j::Int64, g::Float64, b::Float64, ν::Float64)
    inv_dx = 1 / Δx

    # Helper variables.
    # n means next, p means previous.
    h_ij, v1_ij, v2_ij = grid[i, j, :]
    h_in_j, v1_in_j, v2_in_j = grid[i+1, j, :]
    h_ip_j, v1_ip_j, v2_ip_j = grid[i-1, j, :]

    h_i_jn, v1_i_jn, v2_i_jn = grid[i, j+1, :]
    h_i_jp, v1_i_jp, v2_i_jp = grid[i, j-1, :]

    # Derivatives.
    dh_dx = (h_i_jn - h_i_jp) / 2 * inv_dx
    dh_dy = (h_in_j - h_ip_j) / 2 * inv_dx

    dv1_dx = (v1_i_jn - v1_i_jp) / 2 * inv_dx
    dv1_dy = (v1_in_j - v1_ip_j) / 2 * inv_dx

    dv2_dx = (v2_i_jn - v2_i_jp) / 2 * inv_dx
    dv2_dy = (v2_in_j - v2_ip_j) / 2 * inv_dx

    dhv1_dx = (h_i_jn * v1_i_jn - h_i_jp * v1_i_jp) / 2 * inv_dx
    dhv2_dy = (h_in_j * v2_in_j - h_ip_j * v2_ip_j) / 2 * inv_dx

    d2v1_dx2 = (v1_i_jp + v1_i_jn - 2 * v1_ij) * inv_dx^2
    d2v1_dy2 = (v1_ip_j + v1_in_j - 2 * v1_ij) * inv_dx^2

    d2v2_dx2 = (v2_i_jp + v2_i_jn - 2 * v2_ij) * inv_dx^2
    d2v2_dy2 = (v2_ip_j + v2_in_j - 2 * v2_ij) * inv_dx^2

    # RHS formulae.
    rhs = [
        -dhv1_dx - dhv2_dy,
        -g*dh_dx - (v1_ij*dv1_dx + v2_ij*dv1_dy) - b*v1_ij + ν*(d2v1_dx2 + d2v1_dy2),
        -g*dh_dy - (v1_ij*dv2_dx + v2_ij*dv2_dy) - b*v2_ij + ν*(d2v2_dx2 + d2v2_dy2),
    ]

    return rhs
end


function calculate_RHS(grid, Δx::Float64, g::Float64, b::Float64, ν::Float64)
    rhs_grid = zeros(Float64, size(grid))
    # There are ghost cells to implement the boundary conditions.
    # Therefore we omit the first and last rows and columns.
    for j in 2:(size(grid)[2]-1)
        for i in 2:(size(grid)[1]-1)
            rhs_grid[i, j, :] .= calculate_RHS_ij(grid, Δx, i, j, g, b, ν)
        end
    end

    return rhs_grid
end


function explicit_euler_step!(grid, Δt::Float64, Δx::Float64, g::Float64, b::Float64, ν::Float64)
    rhs = calculate_RHS(grid, Δx, g, b, ν)
    grid .+= Δt .* rhs
end


function RK4_step!(grid, Δt::Float64, Δx::Float64, g::Float64, b::Float64, ν::Float64)
    k1 = calculate_RHS(grid                  , Δx, g, b, ν)
    k2 = calculate_RHS(grid .+ Δt .* k1 ./ 2., Δx, g, b, ν)
    k3 = calculate_RHS(grid .+ Δt .* k2 ./ 2., Δx, g, b, ν)
    k4 = calculate_RHS(grid .+ Δt .* k3      , Δx, g, b, ν)
    grid .+= Δt / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
end


function get_north(grid)
    return grid[2, 2:end-1, :]
end
    
function get_south(grid)
    return grid[end-1, 2:end-1, :]
end

function get_west(grid)
    return grid[2:end-1, 2, :]
end

function get_east(grid)
    return grid[2:end-1, end-1, :]
end


function exchange_ghost_cells_serial!(grid)
    # North ghost cells.
    grid[1, 2:end-1, :] .= get_south(grid)
    # South ghost cells.
    grid[end, 2:end-1, :] .= get_north(grid)
    # East ghost cells.
    grid[2:end-1, end, :] .= get_west(grid)
    # West ghost cells.
    grid[2:end-1, 1, :] .= get_east(grid)
end;

##################################################
# Functions for parallel implementation.
##################################################
@everywhere function exchange_ghost_cells_parallel!(grid, neighbours)
    # North ghost cells.
    grid[1, 2:end-1, :] .= remotecall_fetch(get_south, neighbours["top"], grid)
    # South ghost cells.
    grid[end, 2:end-1, :] .= remotecall_fetch(get_north, neighbours["bottom"], grid)
    # East ghost cells.
    grid[2:end-1, end, :] .= remotecall_fetch(get_west, neighbours["right"], grid)
    # West ghost cells.
    grid[2:end-1, 1, :] .= remotecall_fetch(get_east, neighbours["left"], grid)
end;


@everywhere function process_id_to_row(ID, n_cols)
    # The row_ID is in [1, n_rows].
    id_zero = ID - 1
    return Int64(floor(id_zero // n_cols)) + 1
end

@everywhere function process_id_to_col(ID, n_cols)
    # The col_ID is in [1, n_rows].
    id_zero = ID - 1
    return (id_zero % n_cols) + 1
end

@everywhere function shift(ID, period, to_right)
    if to_right
        return (ID % period) + 1
    else
        return ((ID - 2 + period) % period) + 1
    end
end

@everywhere function row_col_to_process_id(row, col, n_cols)
    return (row-1) * n_cols + (col-1) + 1
end


@everywhere function get_neighbours(n_rows, n_cols)
    row_id = process_id_to_row(myid(), n_cols)
    col_id = process_id_to_col(myid(), n_cols)

    right_neighbour = row_col_to_process_id(row_id, shift(col_id, n_cols, true), n_cols)
    left_neighbour = row_col_to_process_id(row_id, shift(col_id, n_cols, false), n_cols)
    top_neighbour = row_col_to_process_id(shift(row_id, n_rows, false), col_id, n_cols)
    bottom_neighbour = row_col_to_process_id(shift(row_id, n_rows, true), col_id, n_cols)
    
    return Dict(
        "right" => right_neighbour,
        "left" => left_neighbour,
        "top" => top_neighbour,
        "bottom" => bottom_neighbour,
    )
end


@everywhere function simulate_parallel(N, L, Δx, Δt, g, b, ν, n_steps, A, σ, n_proc_rows, n_proc_cols)
    ################################
    # Setup process topology
    ################################
    neighbours = get_neighbours(n_proc_rows, n_proc_cols)

    row_id = process_id_to_row(myid(), n_proc_cols) - 1
    col_id = process_id_to_col(myid(), n_proc_cols) - 1

    ################################
    # Distributed initialisation
    ################################
    grid = zeros((
        Int64(ceil(N // n_proc_rows)) + 2,
        Int64(ceil(N // n_proc_cols)) + 2,
        3,
    ))
    # The origin of the grid is the top left corner of this process' subarray.
    grid_origin = [
        -L/2 + L / n_proc_cols * col_id,
        -L/2 + L / n_proc_rows * row_id,
    ]

    init_grid!(grid, grid_origin, Δx, A, σ)

    #################################
    # Distributed evolution
    #################################
    for m in 1:n_steps
        exchange_ghost_cells_parallel!(grid, neighbours)
        RK4_step!(grid, Δt, Δx, g, b, ν)
    end

    return grid
end
