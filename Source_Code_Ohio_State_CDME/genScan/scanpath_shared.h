#pragma once
#include "Layer.h"

/* STRUCTS */
struct SimpleSegment {
	vertex v1, v2;
};

struct BoundingBox {
	vertex tl, tr, bl, br;
};

/* FUNCTIONS */
// Segment
void printSegment(SimpleSegment segment); // Prints a given segment; by personal convention, newlines handled by the caller 
SimpleSegment copySegment(SimpleSegment src); // Deep copy SimpleSegment 
SimpleSegment makeSegment(double x1, double y1, double x2, double y2); // Construct SimpleSegment from x1/y1/x2/y2
inline double Det(double a, double b, double c, double d); // Helper function for lineIntersectsLine
bool getIntersectingVector(SimpleSegment *subsegment, SimpleSegment *segment, BoundingBox *islandBB); // Given a bounding box and a segment, returns true and fills subsegment if there's an overlap, returns false if no overlap
bool lineIntersectsLine(SimpleSegment *s1, SimpleSegment *s2, vertex *intersect); // Tests if lines intersect; if they do, will return true and fill <intersect> with the intersection point; if they don't, returns false

																				  // Vertex
void printVertex(vertex v); // Prints a given vertex; by personal convention, newlines handled by the caller
vertex copyVertex(vertex src); // Deep copy vertex 
bool verticesEqual(vertex v1, vertex v2); // Deep compare vertices 
double distance(vertex v1, vertex v2); // Distance between vertices
void fuzzPoint(vertex *v); // Moves the passed-in vertex slightly; useful to check that you don't have duplicate, identical vectors 

						   // Bounding Box
void printBB(BoundingBox bbox); // Prints a given bbox; by personal convention, newlines handled by the caller
BoundingBox convertBB(vector<vertex> *bbox); // Convert their bbox representation to ours 
bool pointInBB(vertex v1, BoundingBox bb); // Decide if provided point is in bbox

