#include "rayAccelerator.h"
#include "macros.h"
#include <algorithm>

using namespace std;

BVH::BVHNode::BVHNode(void) {}

void BVH::BVHNode::setAABB(AABB& bbox_) { this->bbox = bbox_; }

void BVH::BVHNode::makeLeaf(unsigned int index_, unsigned int n_objs_) {
	this->leaf = true;
	this->index = index_; 
	this->n_objs = n_objs_; 
}

void BVH::BVHNode::makeNode(unsigned int left_index_) {
	this->leaf = false;
	this->index = left_index_; 
			//this->n_objs = n_objs_; 
}


BVH::BVH(void) {}

int BVH::getNumObjects() { return objects.size(); }


void BVH::Build(vector<Object *> &objs) {

		
			BVHNode *root = new BVHNode();

			Vector min = Vector(FLT_MAX, FLT_MAX, FLT_MAX), max = Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX);
			AABB world_bbox = AABB(min, max);

			for (Object* obj : objs) {
				AABB bbox = obj->GetBoundingBox();
				world_bbox.extend(bbox);
				objects.push_back(obj);
			}
			world_bbox.min.x -= EPSILON; world_bbox.min.y -= EPSILON; world_bbox.min.z -= EPSILON;
			world_bbox.max.x += EPSILON; world_bbox.max.y += EPSILON; world_bbox.max.z += EPSILON;
			root->setAABB(world_bbox);
			nodes.push_back(root);
			build_recursive(0, objects.size(), root); // -> root node takes all the 
		}

void BVH::sort(int left_index, int right_index, int max_axis) {
	Comparator comp = Comparator();
	comp.dimension = max_axis;

	std::sort(objects.begin() + left_index, objects.begin() + right_index, comp);
}

void BVH::build_recursive(int left_index, int right_index, BVHNode *node) {
	//PUT YOUR CODE HERE
	if ((right_index - left_index) <= Threshold) {
		// Initiate current node as a leaf with primitives from objects[left_index] to objects[right_index]
		node->makeLeaf(left_index, right_index - left_index);
		return;
	}

	float mid = 0;
	float maxAxis = 0;
	int split_index = 0;

	AABB globalBoundBox = node->getAABB();
	float x = globalBoundBox.max.x - globalBoundBox.min.x;
	float y = globalBoundBox.max.y - globalBoundBox.min.y;
	float z = globalBoundBox.max.z - globalBoundBox.min.z;

	if (x > y && x > z) {
		mid = x;
		maxAxis = 0;
	}

	else if (y > z) {
		mid = y;
		maxAxis = 1;
	}

	else {
		mid = z;
		maxAxis = 2;
	}
	mid = mid / 2;

	
	sort(left_index, right_index, maxAxis);

	for (int i = left_index; i < right_index; i++) {
		if (objects.at(i)->GetBoundingBox().centroid().getAxisValue(i) > mid) {
			split_index = i;
			break;
		}
	}

	if (split_index == left_index || split_index == right_index - 1) {
		int median = (left_index + right_index) / 2;
	
		for (int i = left_index; i < right_index; i++) {
			if (objects.at(i)->GetBoundingBox().centroid().getAxisValue(i) > median) {
				split_index = i;
				break;
			}
		}
	}

	AABB left = AABB(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));
	AABB right = AABB(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));

	// Extend left boundingBox
	for (int i = left_index; i < split_index; i++) {
		AABB box = objects.at(i)->GetBoundingBox();
		left.extend(box);
	}

	// Extend right boundingBox
	for (int i = split_index; i < right_index; i++) {
		AABB box = objects.at(i)->GetBoundingBox();
		right.extend(box);
	}

	// Creating node 
	node->makeNode(left_index);
	BVHNode *left_node = new BVHNode();
	BVHNode *right_node = new BVHNode();
	left_node->setAABB(left);
	right_node->setAABB(right);
	
	nodes.push_back(left_node);
	build_recursive(left_index, split_index, left_node);
	build_recursive(split_index, right_index, right_node);

	//right_index, left_index and split_index refer to the indices in the objects vector
	// do not confuse with left_nodde_index and right_node_index which refer to indices in the nodes vector. 
	// node.index can have a index of objects vector or a index of nodes vector
}

bool BVH::Traverse(Ray& ray, Object** hit_obj, Vector& hit_point) {
			float tmp;
			float tmin = FLT_MAX;  //contains the closest primitive intersection
			bool hit = false;

			BVHNode* currentNode = nodes[0];
			
			//PUT YOUR CODE HERE
			
			return(false);
	}

bool BVH::Traverse(Ray& ray) {  //shadow ray with length
			float tmp;

			double length = ray.direction.length(); //distance between light and intersection point
			ray.direction.normalize();

			//PUT YOUR CODE HERE

			return(false);
	}		
