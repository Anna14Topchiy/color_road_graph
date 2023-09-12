# color_road_graph
A dictionary (road_info_dict) is created to store information about road segments.

A function (are_segments_perpendicular) is defined to check the perpendicularity between two road segments.

The algorithm begins. First, the initial road segment is found, and all connected segments are checked.

If the current segment has no connected segments, a new color is assigned to it.

Then, it is determined which road segments are connected to the current segment, and information about the connected segments is added to the road_info_dict dictionary.

The algorithm continues, and for each segment, its connected segments on the left, right, and in other directions are checked.

Based on the results of the checks, a color is assigned to each segment, and the color_dict dictionary is updated.

A column with colors is added to the GeoDataFrame.

Road segments are displayed on the map using the assigned colors.

The plot is saved as an image with a resolution of 300 dots per inch (dpi), and the image boundaries are compressed to save the entire image to the "solution.png" file.

The image is displayed on the screen.
