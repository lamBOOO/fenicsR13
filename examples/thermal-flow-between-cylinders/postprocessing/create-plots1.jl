using CairoMakie
using CSV
using DataFrames

# Function to load data for a given Kn value
function load_data(Kn)
    return [CSV.File("u_$(Kn)_$i.csv") |> DataFrame for i in 0:3]
end

# Load the data for Kn = 0.05, 0.1, 0.2, 0.4
data_kn_values = Dict(
    0.05 => load_data(0.05),
    0.1  => load_data(0.1),
    0.2  => load_data(0.2),
    0.4  => load_data(0.4)
)

# Colors and line styles for different datasets
line_styles = reverse([:solid, :dash, :dot, :dashdot])
dataset_colors = [:blue, :green, :red, :purple]

# Create a figure with a grid layout: 2 rows x 2 columns
f = Figure(resolution = (800, 600))

# Kn values and their positions in the grid
kn_list = [0.05, 0.1, 0.2, 0.4]
positions = [(1, 1), (1, 2), (2, 1), (2, 2)]

# Iterate over Kn values and their subplot positions
for (Kn_value, pos) in zip(kn_list, positions)
    # Create an axis for each Kn value
    ax = Axis(f[pos...], title = "Kn = $(Kn_value)")
    data_list = data_kn_values[Kn_value]

    # Plot datasets for i = 0:3
    for (i, d) in enumerate(data_list)
        lines!(
            ax,
            d[!, "Points_0"],
            d[!, "u_Magnitude"],
            label = "Mesh $i",
            color = dataset_colors[i],
            linestyle = line_styles[i],
            linewidth = 2
        )
        if i == length(data_dict)
            Legend(f[i, 2], ax)
        end
    end

    # Add legend to the subplot
    # Legend(f[pos[1], pos[2] + 1], ax, title = "Datasets", tellwidth = false)

    # Adjust the subplot layout if needed
    colsize!(f.layout, pos[2], Fixed(300))
    rowsize!(f.layout, pos[1], Fixed(150))
end

save("fig1.png", f)

# Display the figure
f
