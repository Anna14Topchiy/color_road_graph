import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from collections import defaultdict
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


gdf = gpd.read_file('******.shp')

# Create a color dictionary
color_dict = defaultdict(int)

# List of rainbow colors
rainbow_colors = ['#FF0000', '#FF7F00', '#FFFF00', '#00FF00', '#0000FF', '#4B0082', '#8B00FF']

# Create a custom colormap with less saturated colors
num_colors = len(rainbow_colors)
cmap = LinearSegmentedColormap.from_list('custom', rainbow_colors, N=num_colors)

# Criterion for determining matching start and end points of segments
threshold_distance = 1e-5  # Small value for accuracy

# Set equal aspect ratio
plt.gca().set_aspect('equal')

# Create a figure and axes for the plot
fig, ax = plt.subplots(figsize=(10, 10))

# Create a dictionary to store information about roads
road_info_dict = defaultdict(list)

# Function to check perpendicularity of segments
def are_segments_perpendicular(segment1, segment2, threshold_angle):
    if segment1 is not None and segment2 is not None:
        vector1 = np.array(segment1.coords[-1]) - np.array(segment1.coords[0])
        vector2 = np.array(segment2.coords[-1]) - np.array(segment2.coords[0])
        angle_rad = np.arccos(np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2)))
        angle_deg = np.degrees(angle_rad)
        return abs(90 - angle_deg) < threshold_angle
    return False

# Start by finding the initial segment
for index, row in gdf.iterrows():
    current_segment = row['geometry']
    connected_segments = []  # Segments connected to the current one

    for inner_index, inner_row in gdf.iterrows():
        if index != inner_index:
            inner_segment = inner_row['geometry']

            # Check for matching start and end points of the segment with other segments
            if current_segment.distance(inner_segment) < threshold_distance:
                connected_segments.append(inner_segment)

    # If the current segment has no connected segments, consider it as the initial one
    if not connected_segments:
        current_color = len(color_dict)  # Get a new color
        color_dict[index] = current_color
        continue

    # Determine which segments the current segment forms a road with
    road_segments = [current_segment] + connected_segments
    road_info_dict[index] = road_segments  # Add road information to the dictionary

# Continue the algorithm
for index, row in gdf.iterrows():
    current_segment = row['geometry']
    left_connected_segment = None
    right_connected_segment = None
    three_way_connected_segments = []
    two_way_connected_segments = []

    for inner_index, inner_row in gdf.iterrows():
        if index != inner_index:
            inner_segment = inner_row['geometry']

            # Check for matching start and end points of the segment with other segments
            if current_segment.distance(inner_segment) < threshold_distance:
                if are_segments_perpendicular(inner_segment, current_segment, threshold_angle=10):
                    left_connected_segment = inner_segment
                elif are_segments_perpendicular(current_segment, inner_segment, threshold_angle=10):
                    right_connected_segment = inner_segment
                else:
                    two_way_connected_segments.append(inner_segment)

    if left_connected_segment and right_connected_segment:
        # If there are connected segments on the left and right, assign the same color to both
        color_idx = len(color_dict)
        color_dict[index] = color_idx
        color_dict.update({idx: color_idx for idx in gdf.index if
                           idx != index and (gdf.loc[idx, 'geometry'] == left_connected_segment or
                                            gdf.loc[idx, 'geometry'] == right_connected_segment)})
    elif left_connected_segment:
        # If there is only a connected segment on the left, assign it a unique color
        left_color = len(color_dict)
        color_dict[index] = left_color
        color_dict.update({idx: left_color for idx in gdf.index if
                           idx != index and gdf.loc[idx, 'geometry'] == left_connected_segment})
    elif right_connected_segment:
        # If there is only a connected segment on the right, assign it a unique color
        right_color = len(color_dict)
        color_dict[index] = right_color
        color_dict.update({idx: right_color for idx in gdf.index if
                           idx != index and gdf.loc[idx, 'geometry'] == right_connected_segment})
    elif len(two_way_connected_segments) == 2:
        # If there are two connected segments, check their directions (angular coefficients)
        # If they are parallel, assign different colors; otherwise, use the main color
        angle1 = np.arctan2(
            current_segment.coords[1][1] - current_segment.coords[0][1],
            current_segment.coords[1][0] - current_segment.coords[0][0]
        )
        angle2 = np.arctan2(
            two_way_connected_segments[0].coords[1][1] - two_way_connected_segments[0].coords[0][1],
            two_way_connected_segments[0].coords[1][0] - two_way_connected_segments[0].coords[0][0]
        )
        angle3 = np.arctan2(
            two_way_connected_segments[1].coords[1][1] - two_way_connected_segments[1].coords[0][1],
            two_way_connected_segments[1].coords[1][0] - two_way_connected_segments[1].coords[0][0]
        )

        if abs(angle1 - angle2) < np.radians(10) and abs(angle1 - angle3) < np.radians(10):
            left_color = len(color_dict)
            right_color = len(color_dict) + 1

            color_dict[index] = left_color
            color_dict.update({idx: left_color for idx in gdf.index if
                               idx != index and gdf.loc[idx, 'geometry'] == two_way_connected_segments[0]})
            color_dict.update({idx: right_color for idx in gdf.index if
                               idx != index and gdf.loc[idx, 'geometry'] == two_way_connected_segments[1]})
    elif len(two_way_connected_segments) == 1:
        # If there is only one connected segment, assign it a unique color
        connected_color = len(color_dict)
        color_dict[index] = connected_color
        color_dict.update({idx: connected_color for idx in gdf.index if
                           idx != index and gdf.loc[idx, 'geometry'] == two_way_connected_segments[0]})
    elif len(two_way_connected_segments) == 3:
        # If there are three connected segments, use the main color for the current segment
        current_color = color_dict[index]
        color_dict[index] = current_color

# Add a column with colors to the GeoDataFrame
gdf['color'] = gdf.index.map(color_dict)

# Plot segments using the assigned colors
for color in color_dict.values():
    subset = gdf[gdf['color'] == color]
    subset.plot(ax=ax, cmap=cmap, linewidth=1)  # Reduce line thickness for all segments

# Save the plot as an image with specified resolution
plt.savefig('solution.png', dpi=300, bbox_inches='tight')


plt.show()
