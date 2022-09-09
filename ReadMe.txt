
Attached are the files and functions as discussed in my thesis
		---Insert link to page ---

	Code is commented inside files, for information on each function recommend looking there. A brief summary is covered here

	MDS and CSV functions.py
		Accepts the output CSV from the .g files, converts into numpy arrays and, using MDS, plots the distance matrix in Euclidean space.
		Created in python 3.8.5
	
	distance-matrix-generation-functions.g
		File run in GAP systems, required the kbmag package and must be ran in Linux because of this.
		Given a group has a selection of functions to generate a variety of attributes, aimed to estimate the boundary of said group

	Examples and time testing.g
		Gives a brief tutorial on how to use distance-matrix-generation-functions.g
		Along with code used to generate results seen in ---Insert link to pdf ---
		