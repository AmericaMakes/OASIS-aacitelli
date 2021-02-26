// Function definitions and important structs
#include "scanpath_VectorStandardization.h"

// Takes the overall segment, the length of each segment, and the one we want, and figures out both vertices
vertex getScaledPoint(vertex v1, vertex v2, float fractionAlong)
{
	vertex point;
	if (v2.x > v1.x) {
		point.x = v1.x + fractionAlong * (v2.x - v1.x);
	}
	else {
		point.x = v1.x - fractionAlong * (v1.x - v2.x);
	}
	if (v2.y > v1.y) {
		point.y = v1.y + fractionAlong * (v2.y - v1.y);
	}
	else {
		point.y = v1.y - fractionAlong * (v1.y - v2.y);
	}
	return point;
}

// Takes a line segment starting at v1, then pulls it in the direction of v2 for a length of TOO_LONG_CUTOFF. Returns the second endpoint.
vertex extrapolateVector(vertex v1, vertex v2, float TOO_LONG_CUTOFF) {
	// Get angle between points 
	double angle_rads = atan2(v2.y - v1.y, v2.x - v1.x);
	double x_offset = TOO_LONG_CUTOFF * cos(angle_rads);
	double y_offset = TOO_LONG_CUTOFF * sin(angle_rads);

	// Calculate pos of return variable
	vertex v = v1;
	v.x += x_offset;
	v.y += y_offset;
	return v;
}

// Extends the time that short segments take by adding a jump vector after them that makes them take as long as a full-length burn would have 
void lengthenShortVectors(vector<vertex> *list, float TOO_LONG_CUTOFF)
{
	// destList: The list of vectors where the short ones have been lengthened.
	// srcList: The list of vectors we will be lengthening.
	// TOO_LONG_CUTOFF: Anything shorter than this (minus some floating point precision) is considered "too short" and will be lengthened.

	// Iterate through each vector 
	vector<vertex> destList;
	for (int i = 0; i < (*list).size(); i += 2) {

		// Original vector gets added regardless as a mark vector
		vertex v1 = (*list)[i], v2 = (*list)[i + 1];
		destList.push_back(v1);
		destList.push_back(v2);

		// If it's below our goal distance, add a zero-length hatch vector to our extrapolated estimate of where it would be if it were our goal distance
		if (distance(v1, v2) <= TOO_LONG_CUTOFF - .00000001) {
			vertex extrapolatedVertex = extrapolateVector(v1, v2, TOO_LONG_CUTOFF);
			destList.push_back(extrapolatedVertex);
			destList.push_back(extrapolatedVertex);
		}
	}
	*list = destList;
}