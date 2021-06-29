# multiscaleBM

Example of the multiscale BallMapper construction.

Use `multiscale_BM_example.R` to create and save the BM.
MSBM is saved to disk as three csv files:

edges - matrix of a size nr_of_landamrks^2. The value at (i,j) tells at which value of epsilon the the landmark i should be connected with an edge with a landmark j.

distance_of_points_in_order_from_landmarks -- vector of vector (size: number of points times number of landmarks). The vector at i-th
position tells the distances of points from the i-th landmark. All vectors are in ascending ordered.

coloration_in_order_from_landmarks -- as above, the vector tells how the colouration is changing when a new point is added. To get the radius at which change happens, see at the corresponding position at distance_of_points_in_order_from_landmarks.


Use the `multiscale_BM_PLOT.ipynb` notebook to plot the graphs.
Requires Bokeh https://bokeh.org/
