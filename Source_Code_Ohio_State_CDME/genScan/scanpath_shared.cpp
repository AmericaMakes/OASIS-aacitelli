#include "constants.h"
#include "ScanPath.h"
#include "scanpath_shared.h"

/*

SEGMENT

*/
void printSegment(SimpleSegment segment) {
	cout << "(" << segment.v1.x << "," << segment.v1.y << "),(" << segment.v2.x << "," << segment.v2.y << ")";
}

SimpleSegment copySegment(SimpleSegment src) {
	SimpleSegment segment;
	segment.v1.x = src.v1.x;
	segment.v1.y = src.v1.y;
	segment.v2.x = src.v2.x;
	segment.v2.y = src.v2.y;
	return segment;
}

SimpleSegment makeSegment(double x1, double y1, double x2, double y2) {
	SimpleSegment output;
	output.v1.x = x1;
	output.v1.y = y1;
	output.v2.x = x2;
	output.v2.y = y2;
	return output;
}

/*

VERTEX

*/
void printVertex(vertex v) {
	cout << "(" << v.x << "," << v.y << ")";
}

vertex copyVertex(vertex src) {
	vertex v;
	v.x = src.x;
	v.y = src.y;
	return v;
}

bool verticesEqual(vertex v1, vertex v2) {
	return v1.x == v2.x && v1.y == v2.y;
}

double distance(vertex v1, vertex v2) {
	return sqrt(pow(v2.x - v1.x, 2) + pow(v2.y - v1.y, 2));
}

// Used to slightly offset input point to make sure 
void fuzzPoint(vertex *v)
{
	float fuzzAmount = .01 * (((float)rand() / RAND_MAX) * 2 - 1);
	v->x += fuzzAmount;
	v->y += fuzzAmount;
}

/*

BOUNDING BOX

*/
void printBB(BoundingBox bbox) {
	cout << "(";
	printVertex(bbox.tl); cout << ",";
	printVertex(bbox.tr); cout << ",";
	printVertex(bbox.bl); cout << ",";
	printVertex(bbox.br);
	cout << ")";
}

BoundingBox convertBB(vector<vertex> *src) {
	BoundingBox output;
	output.tl.x = (*src)[0].x;
	output.tl.y = (*src)[3].y;
	output.br.x = (*src)[1].x;
	output.br.y = (*src)[2].y;
	output.bl.x = output.tl.x; // Steal these values from ones we already derived, more maintainable
	output.bl.y = output.br.y;
	output.tr.x = output.br.x;
	output.tr.y = output.tl.y;
	return output;
}

bool pointInBB(vertex v1, BoundingBox bb) {
	return v1.x >= bb.tl.x && v1.x <= bb.tr.x && v1.y >= bb.bl.y && v1.y <= bb.tl.y;
}

// Given the overall vector and a given island, returns the subsegment of that vector that is inside the island. 
// THIS IS CALLED WHEN THERE IS NO INTERSECTION. If this is the case, the function will return false.
// Returns true if there was an intersecting vector, false if there wasn't.
bool getIntersectingVector(SimpleSegment *subsegment, SimpleSegment *segment, BoundingBox *islandBB)
{
	// Default is returning two vertices with null values.
	subsegment->v1.x = NULL;
	subsegment->v1.y = NULL;
	subsegment->v2.x = NULL;
	subsegment->v2.y = NULL;

	// Construct segments for each side of our island (for comparison purposes)
	SimpleSegment top, left, right, bottom;
	top.v1 = islandBB->tl; top.v2 = islandBB->tr;
	left.v1 = islandBB->tl; left.v2 = islandBB->bl;
	right.v1 = islandBB->tr; right.v2 = islandBB->br;
	bottom.v1 = islandBB->bl; bottom.v2 = islandBB->br;

	// Edge Case: Colinear to one of the edges. The function implementations I am using are a little undefined here so I handle them manually to be sure.
	// This will likely never happen due to floating point precision, but gotta be ready incase it does! 

	// IF segment is straight horizontal AND segment is straight horizontal to one of the vertical portions
	// Note that we don't consider anything right colinear, only left, to avoid duplicates entries across islands.
	if (segment->v1.y == segment->v2.y && segment->v1.y == islandBB->bl.y) {

		// If there is no overlap, get out.
		if ((segment->v1.x <= islandBB->tl.x && segment->v2.x <= islandBB->tl.x) ||
			segment->v1.x >= islandBB->tl.x && segment->v2.x >= islandBB->tl.x) {
			return false;
		}

		vertex leftPoint = segment->v1.x <= segment->v2.x ? segment->v1 : segment->v2;
		vertex rightPoint = segment->v1.x <= segment->v2.x ? segment->v2 : segment->v1;

		// Figure out left endpoint 
		double left;
		if (leftPoint.x >= islandBB->tl.x) {
			left = leftPoint.x;
		}
		else {
			left = islandBB->tl.x;
		}

		// Determine right endpoint
		double right;
		if (rightPoint.x <= islandBB->tr.x) {
			right = rightPoint.x;
		}
		else {
			right = islandBB->tr.x;
		}

		// Subsegment is the y we already established, then left to right boundaries
		*subsegment = makeSegment(left, segment->v1.y, right, segment->v1.y);
		return true;
	}

	// Specifically want to exclude top ones as previously mentioned (to avoid duplicates)
	if (segment->v1.y == segment->v2.y && segment->v1.y == islandBB->tl.y) {
		return false;
	}

	// Vertical Cases
	// Note that we don't consider anything top colinear, only bottom, to avoid duplicate entries across islands.
	if (segment->v1.x == segment->v2.x && segment->v1.x == islandBB->tl.x) {

		// If there is no overlap, get out.
		if ((segment->v1.y <= islandBB->bl.y && segment->v2.y <= islandBB->bl.y) ||
			segment->v1.y >= islandBB->tl.y && segment->v2.y >= islandBB->tl.y) {
			return false;
		}

		vertex bottomPoint = segment->v1.y <= segment->v2.y ? segment->v1 : segment->v2;
		vertex topPoint = segment->v1.y <= segment->v2.y ? segment->v2 : segment->v1;

		// Figure out left endpoint 
		double bottom;
		if (bottomPoint.y >= islandBB->bl.y) {
			bottom = bottomPoint.y;
		}
		else {
			bottom = islandBB->bl.y;
		}

		// Determine right endpoint
		double top;
		if (topPoint.y <= islandBB->tl.y) {
			top = topPoint.y;
		}
		else {
			top = islandBB->tl.y;
		}

		// Subsegment is the y we already established, then left to right boundaries
		*subsegment = makeSegment(segment->v1.x, bottom, segment->v1.x, top);
		return true;
	}

	// Specifically want to exclude right ones as previously mentioned (to avoid duplicates)
	if (segment->v1.x == segment->v2.x && segment->v1.x == islandBB->tr.x) {
		return false;
	}

	// Check each side. If there is an intersection we know the intersection point will be filled, so we can then check that it's actually within the segment.
	// Again, note that lineIntersectsLine just returns if they *ever* intersect, and if they do it returns a point. If we know they intersect we have to compound it with checking that it's within our bounds.
	vertex topIntersection, rightIntersection, bottomIntersection, leftIntersection;
	bool intersectsTop = lineIntersectsLine(&top, segment, &topIntersection) && // They aren't parallel
		(topIntersection.x >= islandBB->tl.x && topIntersection.x <= islandBB->tr.x) && // Their intersection point is on the actual island boundary itself 
		((segment->v1.y <= islandBB->tl.y && segment->v2.y >= islandBB->tl.y) || // This and last condition ensure segment is actually split across, not just intersecting if we extend it
		(segment->v1.y >= islandBB->tl.y && segment->v2.y <= islandBB->tl.y));
	bool intersectsRight = lineIntersectsLine(&right, segment, &rightIntersection) &&
		(rightIntersection.y >= islandBB->bl.y && rightIntersection.y <= islandBB->tl.y) &&
		((segment->v1.x <= islandBB->tr.x && segment->v2.x >= islandBB->tr.x) ||
		(segment->v1.x >= islandBB->tr.x && segment->v2.x <= islandBB->tr.x));
	bool intersectsBottom = lineIntersectsLine(&bottom, segment, &bottomIntersection) &&
		(bottomIntersection.x >= islandBB->tl.x && bottomIntersection.x <= islandBB->tr.x) &&
		((segment->v1.y <= islandBB->bl.y && segment->v2.y >= islandBB->bl.y) ||
		(segment->v1.y >= islandBB->bl.y && segment->v2.y <= islandBB->bl.y));
	bool intersectsLeft = lineIntersectsLine(&left, segment, &leftIntersection) &&
		(leftIntersection.y >= islandBB->bl.y && leftIntersection.y <= islandBB->tl.y) &&
		((segment->v1.x <= islandBB->tl.x && segment->v2.x >= islandBB->tl.x) ||
		(segment->v1.x >= islandBB->tl.x && segment->v2.x <= islandBB->tl.x));

	int numIntersections = 0;
	if (intersectsTop) numIntersections++;
	if (intersectsRight) numIntersections++;
	if (intersectsBottom) numIntersections++;
	if (intersectsLeft) numIntersections++;

	// If no intersections
	if (numIntersections == 0) {

		// It's either both are in or neither is in at this point. If the first is out, whole segment is out.
		if (segment->v1.x < islandBB->tl.x || segment->v1.x > islandBB->tr.x ||
			segment->v1.y < islandBB->bl.y || segment->v1.y > islandBB->tl.y) {
			return false;
		}

		// Otherwise, whole segment is in.
		*subsegment = copySegment(*segment); // Normally I'd just switch the pointer but we want to copy it all over to avoid aliasing 
	}

	// If one intersection, the line is the one intersection point we know plus whichever endpoint is within the rectangle.
	else if (numIntersections == 1) {

		// Figure out which intersection is inside the rectangle
		if (segment->v1.x >= islandBB->tl.x && segment->v1.x <= islandBB->tr.x &&
			segment->v1.y >= islandBB->bl.y && segment->v1.y <= islandBB->tl.y) {
			subsegment->v1 = copyVertex(segment->v1);
		}
		else {
			subsegment->v1 = copyVertex(segment->v2);
		}

		// Other half of the vector is from our original intersection calculations
		if (intersectsTop) {
			subsegment->v2 = copyVertex(topIntersection);
		}
		else if (intersectsRight) {
			subsegment->v2 = copyVertex(rightIntersection);
		}
		else if (intersectsBottom) {
			subsegment->v2 = copyVertex(bottomIntersection);
		}
		else if (intersectsLeft) {
			subsegment->v2 = copyVertex(leftIntersection);
		}
	}

	// If two intersections, it must intersect two sides of our rectangle. We add any of the two vertices that aren't duplicates.
	// Also covers 3 and 4 intersection cases because that just means a vertex double counts. We just make sure we pick two vertices of our set that aren't the same.
	else if (numIntersections >= 2 && numIntersections <= 4) {

		// Edge Case: One point on island corner and other point outside of island in a direction that has no overlap w/ island
		if (numIntersections == 2) {

			// Top-Right: Error out in case that intersection is exactly on and other point is outside island
			if (intersectsTop && intersectsRight &&
				((verticesEqual(segment->v1, islandBB->tr) && !pointInBB(segment->v2, *islandBB)) ||
				(verticesEqual(segment->v2, islandBB->tr) && !pointInBB(segment->v1, *islandBB)))) {
				return false;
			}

			// Check top-left
			if (intersectsTop && intersectsLeft &&
				((verticesEqual(segment->v1, islandBB->tl) && !pointInBB(segment->v2, *islandBB)) ||
				(verticesEqual(segment->v2, islandBB->tl) && !pointInBB(segment->v1, *islandBB)))) {
				return false;
			}

			// Check bottom-right
			if (intersectsBottom && intersectsRight &&
				((verticesEqual(segment->v1, islandBB->br) && !pointInBB(segment->v2, *islandBB)) ||
				(verticesEqual(segment->v2, islandBB->br) && !pointInBB(segment->v1, *islandBB)))) {
				return false;
			}

			// Check bottom-left 
			if (intersectsBottom && intersectsLeft &&
				((verticesEqual(segment->v1, islandBB->bl) && !pointInBB(segment->v2, *islandBB)) ||
				(verticesEqual(segment->v2, islandBB->bl) && !pointInBB(segment->v1, *islandBB)))) {
				return false;
			}
		}

		bool firstFilled = false, secondFilled = false;

		// Handle top intersection 
		if (intersectsTop) {
			subsegment->v1 = copyVertex(topIntersection);
			firstFilled = true;
		}

		// Handle bottom intersection 
		if (intersectsBottom) {

			// See if first AND check that it's not a duplicate. 
			// Note that short-circuiting results in this never accessing something that hasn't been assigned its value in our system yet.
			if (firstFilled && !verticesEqual(subsegment->v1, bottomIntersection)) {
				subsegment->v2 = copyVertex(bottomIntersection);
				secondFilled = true;
			}

			else if (!firstFilled) {
				subsegment->v1 = copyVertex(bottomIntersection);
				firstFilled = true;
			}
		}

		// Handle right intersection
		if (!secondFilled && intersectsRight) {

			// See if second AND check that it's not a duplicate
			if (firstFilled && !verticesEqual(subsegment->v1, rightIntersection)) {
				subsegment->v2 = copyVertex(rightIntersection);
				secondFilled = true;
			}

			// Otherwise, see if first
			else if (!firstFilled)
			{
				subsegment->v1 = copyVertex(rightIntersection);
				firstFilled = true;
			}
		}

		// Handle left intersection
		if (!secondFilled && intersectsLeft) {

			// See if second AND check that it's not a duplicate
			if (firstFilled && !verticesEqual(subsegment->v1, leftIntersection)) {
				subsegment->v2 = copyVertex(leftIntersection);
				secondFilled = true;
			}

			// Otherwise, see if first
			else if (!firstFilled)
			{
				subsegment->v1 = copyVertex(leftIntersection);
				firstFilled = true;
			}
		}

		// Error case that provides more debug information
		if (!secondFilled) {
			cout << "getIntersectingVector: Got through 2/3/4 intersection code without two intersections!" << endl;
			cout << "numIntersections: " << numIntersections << endl;
			cout << "topIntersection: "; printVertex(topIntersection); cout << endl;
			cout << "bottomIntersection: "; printVertex(bottomIntersection); cout << endl;
			cout << "leftIntersection: "; printVertex(leftIntersection); cout << endl;
			cout << "rightIntersection: "; printVertex(rightIntersection); cout << endl;
			cout << "segment: "; printSegment(*segment); cout << endl;
			cout << "bb: "; printBB(*islandBB); cout << endl;
			return false;
		}
	}

	// Error case that provides more debug information
	else {
		cout << "getIntersectingVector: numIntersections was " << numIntersections << " instead of 0/1/2/3!" << endl;
		return false;
	}

	// Assume we found one; We return false elsewhere if we hit any of the error conditions
	return true;
}

inline double Det(double a, double b, double c, double d)
{
	return a*d - b*c;
}

/// Calculate intersection of two lines.
///\return true if found, false if not found or error
bool lineIntersectsLine(SimpleSegment *s1, SimpleSegment *s2, vertex *intersect)
{
	// Found this implementation online; easier to just translate my system of segments here than replace every one and risk screwing it up
	double x1 = s1->v1.x, y1 = s1->v1.y;
	double x2 = s1->v2.x, y2 = s1->v2.y;
	double x3 = s2->v1.x, y3 = s2->v1.y;
	double x4 = s2->v2.x, y4 = s2->v2.y;


	double detL1 = Det(x1, y1, x2, y2);
	double detL2 = Det(x3, y3, x4, y4);
	double x1mx2 = x1 - x2;
	double x3mx4 = x3 - x4;
	double y1my2 = y1 - y2;
	double y3my4 = y3 - y4;

	double denom = Det(x1mx2, y1my2, x3mx4, y3my4);
	if (denom == 0.0) // Lines don't seem to cross
	{
		intersect->x = NAN;
		intersect->y = NAN;
		return false;
	}

	double xnom = Det(detL1, x1mx2, detL2, x3mx4);
	double ynom = Det(detL1, y1my2, detL2, y3my4);
	intersect->x = xnom / denom;
	intersect->y = ynom / denom;
	if (!isfinite(intersect->x) || !isfinite(intersect->y)) // Probably a numerical issue
		return false;

	return true; //All OK
}