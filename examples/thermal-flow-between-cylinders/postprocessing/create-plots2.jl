using CairoMakie
using CSV
using DataFrames

# Function to load data for a given Kn value
function load_data(Kn)
    return [CSV.File("u_$(Kn)_$i.csv") |> DataFrame for i in 0:3]
end

# Load the data for Kn = 0.05, 0.1, 0.2, 0.4
data_0_05 = load_data(0.05)
data_0_1  = load_data(0.1)
data_0_2  = load_data(0.2)
data_0_4  = load_data(0.4)

# Create a figure
f = Figure()

# Create an axis
ax1 = Axis(f[1, 1], title = "Velocity Magnitude vs. Position w.r.t. mesh sizes", aspect = 1)

# Colors for different Kn values
kn_colors = Dict(
    0.05 => :blue,
    0.1  => :green,
    0.2  => :red,
    0.4  => :purple
)

# Plotting function
function plot_data(ax, data_list, Kn_value, color)
    for (i, d) in enumerate(data_list)
        lines!(
            ax,
            d[!, "Points_0"],
            d[!, "u_Magnitude"],
            label = "Kn = $(Kn_value), Mesh $i",
            color = color,
            linestyle = [:solid, :dash, :dot, :dashdot][i]
        )
    end
end

# Plot data for each Kn value
plot_data(ax1, data_0_05, 0.05, kn_colors[0.05])
plot_data(ax1, data_0_1,  0.1,  kn_colors[0.1])
plot_data(ax1, data_0_2,  0.2,  kn_colors[0.2])
plot_data(ax1, data_0_4,  0.4,  kn_colors[0.4])

# Create a legend
legend = Legend(f[1, 2], ax1)

# Adjust the layout
colsize!(f.layout, 1, Fixed(300))  # Adjust as needed

save("fig2.png", f)

# Display the figure
f

