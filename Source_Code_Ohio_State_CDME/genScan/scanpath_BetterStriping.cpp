// Function definitions and important structs
#include "scanpath_BetterStriping.h"

// Gives us access to M_PI and trig functions
#define _USE_MATH_DEFINES
#include <math.h>

// Represent line as one (x, y) point and a slope
// Note that the range [89deg, 91deg] is collapsed to 90deg to avoid integer overflows, and the slope is set to INT_MAX, which is a special case in many places
struct Line {
	double x, y, m; 
};

// In order to support "ranges" that are at angles other than 0/90/180/270, we instead represent areas as two parallel lines, or an "alley"
struct Alley {
	Line l1;
	Line l2;
};

// Given a point and an angle, gives its y-intercept
// 99% SURE THIS WORKS 
double getYInt(double x, double y, double angle_rads) {
	return y - tan(angle_rads) * x; // Just a reformulation of slope-intercept form; has to convert angle to slope via tan
}

// Given a point and an angle, give its x-intercept
// 99% SURE THIS WORKS 
double getXInt(double x, double y, double angle_rads) {

	// Definitely a way to do this without finding y intercept first, but it's been like eight years since I worked with all this slope stuff so cut me some slack 
	double yInt = getYInt(x, y, angle_rads);

	// 0 = mx + b
	// -b = mx
	// x = -b/m
	return -1 * yInt / tan(angle_rads);
}

// Gives the y-value of a given Line at the given x value 
// 99% SURE THIS WORKS 
double valueAtX(double x, Line l) {
	return l.y + l.m * (x - l.x); // Essentially slope intercept but with the provided point as the origin
}

// Convert degrees to rads
double rads(double degs) {
	return degs / 180 * M_PI;
}

/* Generates the alleys used to split the vectors along.
	My way of doing this is essentially to generate vertical lines stripeWidth apart, tilt them by hatchAngle. 
	Note that, if we just do that, though, the actual distance between them will be less by some small number you can get with trig.
	Therefore, I initially generate the lines stripeWidth * <that number> apart so that, once we rotate them, they are stripeWidth apart.
	Also note that this works differently depending which direction the hatches are (x/y). It's the same problem solved by hatch(), so we just pass in the result from that. */
void generateAlleys(vector<Alley> *alleys, string primaryHatchDir, double alleyAngle, double stripeWidth, BoundingBox *bbox) {

	// If hatches are vertical, we generate horizontal alleys by spacing start points in the vertical direction 
	if (primaryHatchDir.compare("x") == 0) {

		// 1. Determine begin/end points such that the tilting doesn't introduce missed areas
		// If the angle is downward, we would need to start upward by (bbox width / tan(hatchAngle)) so there isn't uncovered area in the top-right of the part
		// Rather than do some sketchy trig, this is low enough performance gain that we can just make a conservative estimate to pad on each side
		// Via trig this should just be tan(45deg) which is just 1, so a full width onto it, but to be safe I multiply it by 2
		// Can revert if it becomes a legitimate performance penalty once we iterate bigger
		double start = bbox->bl.y - (bbox->tr.x - bbox->tl.x);
		/*
		double start = bbox->bl.y;
		if (alleyAngle >= 90 && alleyAngle <= 180) { // Need this case because, for example, 30deg and 210deg are visually the same line but they'll mess up the math if I don't account for them 
			double angle = abs(180 - alleyAngle); // Gets angle from the 180deg line, which is simpler conceptually for me to do math with 
			start -= (bbox->tr.x - bbox->tl.x) * abs(sin(rads(angle))); // Divide by 2 is because we are
		}
		*/

		// If the angle is more upward, we would need to correct downward by (bbox width  * tan(hatchAngle)) so there isn't uncovered area in the bottom-right of the part 
		double end = bbox->tl.y + (bbox->tr.x - bbox->tl.x);
		/*
		double end = bbox->tl.y;
		if (alleyAngle < 90 && alleyAngle >= 0) {
			end += (bbox->tr.x - bbox->tl.x) * abs(sin(rads(alleyAngle)));
		}
		*/

		// 2. Calculate spacing such that tilting brings it to our desired stripe width, not less than it 
		double spacing = stripeWidth / abs(cos(rads(alleyAngle)));
		// cout << "given hatchAngle " << hatchAngle << " and stripeWidth " << stripeWidth << ", spacing = " << spacing << endl;

		// 3. Actually generate the lines 
		for (double i = start; i < end; i += spacing) { // Space from <bbox.minY> to <bbox.maxY> (plus one more region in order to fully overlap the original bbox) by <stripeWidth * magicNumber> space 
			
			Alley alley;
			
			// Derive slope from angle (we know it's [0, 45] U [135, 180] so we don't have to check for 90deg)
			alley.l1.m = tan(rads(alleyAngle));
			alley.l2.m = tan(rads(alleyAngle));

			// We store these with their y-intercepts; if we have a line of slope zero, storing it as its x-intercept gets really messy
			alley.l1.x = alley.l2.x = bbox->bl.x;
			alley.l1.y = i;
			alley.l2.y = i + spacing; 
			
			(*alleys).push_back(alley);
		}
	}

	// If hatches are horizontal, we generate vertical alleys by spacing start points in the horizontal direction
	else {
		
		// 1. Determine begin/end points such that the tilting doesn't introduce missed areas 
		// If the angle is rightward, we need to start leftward by () so that we don't miss space in top-left of part  
		double start = bbox->tl.x - (bbox->tl.y - bbox->bl.y);
		/*
		double start = bbox->tl.x;
		if (alleyAngle >= 0 && alleyAngle <= 90) {
			start -= (bbox->tl.y - bbox->bl.y) * abs(cos(rads(alleyAngle)));
		}
		*/

		// If the angle is leftward, we need to start rightward by () so that we don't miss space in top-right of part 
		double end = bbox->tr.x + (bbox->tl.y - bbox->bl.y);
		/*
		double end = bbox->tr.x;
		if (alleyAngle > 90 && alleyAngle < 180) {
			double angle = alleyAngle - 90;
			end += (bbox->tl.y - bbox->bl.y) * abs(cos(rads(angle)));
		}
		*/

		// 2. Calculate spacing such that tilting brings it to our desired stripe width, not less than it 
		double spacing = stripeWidth / abs(sin(rads(alleyAngle)));

		// 3. Actually generate the lines
		for (double i = start; i < end; i += spacing) {
			Alley alley;

			// Derive slope from hatchAngle
			// We know alleyAngle is [0, 89] U [90] U [91, 180], so this is the only edge case / integer overflow check we have to do
			if (alleyAngle != 90) {
				alley.l1.m = tan(rads(alleyAngle));
				alley.l2.m = tan(rads(alleyAngle));
				alley.l1.x = i;
				alley.l1.y = bbox->bl.y;
				alley.l2.x = i + spacing;
				alley.l2.y = bbox->bl.y;
			}
			else {
				alley.l1.m = INT_MAX; // (tangent undefined at 90deg)
				alley.l2.m = INT_MAX;
				alley.l1.x = i;
				alley.l2.x = i + spacing;

				// y-coord is redundant on these; unfortunately C++ won't error if we try to access it (it'll just interpret whatever garbage is in memory as a number) so we just set it to 0
				alley.l1.y = alley.l2.y = 0;
			}

			(*alleys).push_back(alley);
		}
	}
}

// Source: https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line (essentially does a cross product)
// Negative if the point is to the left of the line (or above if the line is 0deg)
// Zero if the point is on the line 
// Positive if the point is to the left of the line (or below if the line is 0deg)
int sideOfLine(vertex a, vertex b, vertex c) {

	// Order points lexicographically so that return of this function doesn't change based on which order you provide a,b in
	if (a.y > b.y) {
		vertex swap;
		swap = a;
		a = b;
		b = swap;
	}

	return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// Reverses in groups of two to keep vector directions the same
// I manually examined the output of this and it seems to work as intended 
void reverseVertexVector(vector<vertex> *v) {
	vector<vertex> swap;
	for (int i = (*v).size() - 2; i >= 0; i -= 2) {
		swap.push_back((*v)[i]);
		swap.push_back((*v)[i + 1]);
	}
	(*v) = swap;
}

// Converts a Line to a SimpleSegment (I end up using this in various places, might as well turn it into a function)
SimpleSegment makeSegmentFromLine(Line l) {
	SimpleSegment output;
	if (l.m != INT_MAX) {
		output.v1.x = -1000;
		output.v1.y = valueAtX(output.v1.x, l);
		output.v2.x = 1000;
		output.v2.y = valueAtX(output.v2.x, l);
	}
	else {
		output.v1.x = output.v2.x = l.x;
		output.v1.y = -1000;
		output.v2.y = 1000;
	}
	return output;
}

// lineIntersectsLine returns true in any circumstance they aren't parallel;
// This function, when provided that intersection and that segment, returns whether that intersection is actually within the bounds of the segment 
bool intersectionWithinBounds(vertex intersection, SimpleSegment segment) {
	return intersection.x >= min(segment.v1.x, segment.v2.x) && intersection.x <= max(segment.v1.x, segment.v2.x) &&
		intersection.y >= min(segment.v1.y, segment.v2.y) && intersection.y <= max(segment.v1.y, segment.v2.y);
}

// If an overlap exists between the specified segment and alley, fills <overlap> and returns true.
// If an overlap does not exist between the specified segment and alley, returns false.
bool getSegmentAlleyIntersection(SimpleSegment *overlap, SimpleSegment segment, Alley alley, string primaryHatchDir) {

	// Use this a ton, might as well calculate it once 
	SimpleSegment leftLineSegment = makeSegmentFromLine(alley.l1);
	SimpleSegment rightLineSegment = makeSegmentFromLine(alley.l2);

	// Special case for vertical alleys 
	if (alley.l1.m == INT_MAX) {

		// Case 1: Both are left or both are right, meaning we return false and don't fill overlap 
		if (segment.v1.x < alley.l1.x && segment.v2.x < alley.l1.x ||
			segment.v1.x > alley.l2.x && segment.v2.x > alley.l2.x) {
			return false;
		}

		// Case 2: One is left and one is inside, meaning our subset is from the left intersection to the point inside 
		if (segment.v1.x < alley.l1.x && segment.v2.x > alley.l1.x && segment.v2.x < alley.l2.x ||
			segment.v2.x < alley.l1.x && segment.v1.x > alley.l1.x && segment.v1.x < alley.l2.x) {
			vertex intersect;
			lineIntersectsLine(&leftLineSegment, &segment, &intersect);
			(*overlap).v1 = intersect;
			(*overlap).v2 = segment.v1.x > segment.v2.x ? segment.v1 : segment.v2; // Pick rightmost one 
			return true;
		}

		// Case 3: One is right and one is inside, meaning our subset is from the right intersection to the point inside 
		if (segment.v1.x > alley.l1.x && segment.v1.x < alley.l2.x && segment.v2.x > alley.l2.x ||
			segment.v2.x > alley.l1.x && segment.v2.x < alley.l2.x && segment.v1.x > alley.l2.x) {
			vertex intersect;
			lineIntersectsLine(&rightLineSegment, &segment, &intersect);
			(*overlap).v1 = segment.v1.x < segment.v2.x ? segment.v1 : segment.v2; // Pick leftmost one 
			(*overlap).v2 = intersect;
			return true;
		}

		// Case 4: Both are inside, meaning we just return our original vector 
		if (segment.v1.x > alley.l1.x && segment.v1.x < alley.l2.x &&
			segment.v2.x > alley.l1.x && segment.v2.x < alley.l2.x) {
			*overlap = segment;
			return true;
		}

		// Case 5: They are split across, meaning we set equal to both intersection points 
		if (segment.v1.x < alley.l1.x && segment.v2.x > alley.l2.x ||
			segment.v1.x > alley.l2.x && segment.v2.x < alley.l1.x) {
			lineIntersectsLine(&leftLineSegment, &segment, &((*overlap).v1));
			lineIntersectsLine(&rightLineSegment, &segment, &((*overlap).v2));
			return true;
		}

		// Logically should never get here
		else {
			cout << "getSegmentAlleyIntersection hit an unrecognized vertical-line case!" << endl;
			return false;
		}
	}

	// See if segment intersects 0, 1, or 2 times
	// NOTE: lineIntersectsLine returns true in any circumstance that the lines aren't parallel, which is why we need the extra checks 
	vertex leftIntersection, rightIntersection;
	bool leftIntersects = lineIntersectsLine(&leftLineSegment, &segment, &leftIntersection) && intersectionWithinBounds(leftIntersection, segment);
	bool rightIntersects = lineIntersectsLine(&rightLineSegment, &segment, &rightIntersection) && intersectionWithinBounds(rightIntersection, segment);

	// If 0 intersections, it's either completely inside or completely outside. We test one of the endpoints; if inside, whole subsegment. If outside, no subsegment. 
	if (!leftIntersects && !rightIntersects) {
		*overlap = segment; 
		if (sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v1) > 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v1) < 0 ||
			sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v1) < 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v1) > 0) {
			return true;
		}
		return false;
	}

	else if (leftIntersects != rightIntersects) {
		if (leftIntersects) {
			(*overlap).v2 = leftIntersection;
			if (sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v1) > 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v1) < 0 ||
				sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v1) < 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v1) > 0) {
				(*overlap).v1 = segment.v1;
			}
			else if (sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v2) > 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v2) < 0 ||
					 sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v2) < 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v2) > 0) {
				(*overlap).v1 = segment.v2;
			}
			else { // One of the endpoints is exactly on an alley (should never happen b/c floating point, but here just incase
				return false;
			}
			return true;
		}

		else {
			(*overlap).v1 = rightIntersection;
			if (sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v1) > 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v1) < 0 ||
				sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v1) < 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v1) > 0) {
				(*overlap).v2 = segment.v1;
			}
			else if (sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v2) > 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v2) < 0 ||
					 sideOfLine(leftLineSegment.v1, leftLineSegment.v2, segment.v2) < 0 && sideOfLine(rightLineSegment.v1, rightLineSegment.v2, segment.v2) > 0) {
				(*overlap).v2 = segment.v2;
			}
			else { // One of the endpoints is exactly on an alley (should never happen b/c floating point, but here just incase
				return false;
			}
			return true;
		}
	}

	// If 2 intersections, subsegment is from one intersection to the other 
	else {
		(*overlap).v1 = leftIntersection;
		(*overlap).v2 = rightIntersection;
		return true;
	}

	// Code returns in every circumstance, don't need another return 
}
 
/* Essentially does the following:
	1. Determine our striping ranges
	2. Foreach striping range, iterate through every existing vector, strip out the section that is within the striping range's bounds, then add it to a separate output queue.
		This is how we reorder them in the order they show up in the striping alleys.
*/
void betterStriping(vector<vertex> *isList, vector<vertex> *bbox, string primaryHatchDir, double hatchAngle, double stripeWidth) {

	// 1. Convert their bbox representation to ours
	BoundingBox boundingBox = convertBB(bbox);

	// 2. Generate our alleys; all we need for this is hatchAngle and stripeWidth
	vector<Alley> alleys;
	double alleyAngle = fmod(hatchAngle + 90, 180);
	alleyAngle = alleyAngle >= 88 && alleyAngle <= 92 ? 90 : alleyAngle; // Round values close to 90 to 90 to avoid integer overflow circumstances
	generateAlleys(&alleys, primaryHatchDir, alleyAngle, stripeWidth, &boundingBox);

	// NOTE: If alleys are vertical, alleys *should* already be sorted by x-coordinate. We want to progress through alleys left to right for gas flow direction reasons.

	// 3. Iterate through alleys, taking out the subset of each vector that falls in each alley
	vector<vertex> output, outputSubset;
	for (int i = 0; i < alleys.size(); i++) {
		for (int j = 0; j < (*isList).size(); j += 2) {
			SimpleSegment segment;
			segment.v1 = (*isList)[j];
			segment.v2 = (*isList)[j + 1];

			// Check for overlap and add to output list if overlap exists
			SimpleSegment overlap;
			if (getSegmentAlleyIntersection(&overlap, segment, alleys[i], primaryHatchDir)) {
				outputSubset.push_back(overlap.v1);
				outputSubset.push_back(overlap.v2);
			}
			// If line doesn't overlap with given alley, we don't add anything to the output list 
		}
		
		// We want to sort this alley's subsets from min-x to max-x (REGARDLESS OF HORIZONTAL OR VERTICAL) due to gas flow direction 
		if (outputSubset.size() > 0) {
			// We already know they're sorted one way, so we just check if the way they are sorted is correct and reverse it if they aren't
			if (outputSubset[0].x > outputSubset[outputSubset.size() - 1].x) {
				reverseVertexVector(&outputSubset);
			}
		}

		// After we've given the chance for the code to sort by x-coordinate, move them to our final output list 
		for (int j = 0; j < outputSubset.size(); j++) {
			output.push_back(outputSubset[j]);
		}
		outputSubset.clear();
	}

	// Copy our output back into isList
	// There exists a better way to do this with pointers but this is much easier to program 
	*isList = output;

	// Fuzz all the points prior to our alley overlay (helps identify duplicates)
	/*
	for (int i = 0; i < (*isList).size(); i++) {
		fuzzPoint(&((*isList)[i]));
	}
	*/

	// Draw Valleys (PURELY FOR DEBUGGING; thus, I won't worry about making it perfect for anything more than my test cases)
	/*
	for (int i = 0; i < alleys.size(); i++) {
		SimpleSegment leftLine = makeSegmentFromLine(alleys[i].l1);
		(*isList).push_back(leftLine.v1);
		(*isList).push_back(leftLine.v2);
	}
	*/
}