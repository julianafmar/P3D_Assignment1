#include "rayAccelerator.h"
#include "macros.h"
#include <algorithm>

using namespace std;

int maxAxis = 0;

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


void BVH::Build(vector<Object*>& objs) {


	BVHNode* root = new BVHNode();

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

void BVH::sort(int left_index, int right_index) {
	Comparator comp = Comparator();
	comp.dimension = maxAxis;

	std::sort(objects.begin() + left_index, objects.begin() + right_index, comp);
}

void BVH::build_recursive(int left_index, int right_index, BVHNode* node) {
	if ((right_index - left_index) <= Threshold) {
		// Initiate current node as a leaf with primitives from objects[left_index] to objects[right_index]
		node->makeLeaf(left_index, right_index - left_index);
		return;
	}

	float mid = 0;
	int split_index = left_index;

	AABB globalBoundBox = node->getAABB();
	float x = globalBoundBox.max.x - globalBoundBox.min.x;
	float y = globalBoundBox.max.y - globalBoundBox.min.y;
	float z = globalBoundBox.max.z - globalBoundBox.min.z;

	// Find the axis with the maximum length

	if (x >= y) {
		mid = x;
		maxAxis = 0;
	}

	else {
		mid = y;
		maxAxis = 1;
	}

	if (z > mid) {
		mid = z;
		maxAxis = 2;
	}
	mid = mid / 2;
	mid += globalBoundBox.min.getAxisValue(maxAxis);


	sort(left_index, right_index);

	for (int i = left_index; i < right_index; i++) {
		if (objects.at(i)->getCentroid().getAxisValue(maxAxis) > mid) {
			split_index = i;
			break;
		}
	}

	if (split_index == left_index || split_index == right_index - 1) {
		split_index = (left_index + right_index) / 2;

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
	BVHNode* left_node = new BVHNode();
	BVHNode* right_node = new BVHNode();

	node->makeNode(nodes.size());
	left_node->setAABB(left);
	right_node->setAABB(right);

	nodes.push_back(left_node);
	nodes.push_back(right_node);
	build_recursive(left_index, split_index, left_node);
	build_recursive(split_index, right_index, right_node);

	//right_index, left_index and split_index refer to the indices in the objects vector
	// do not confuse with left_nodde_index and right_node_index which refer to indices in the nodes vector. 
	// node.index can have a index of objects vector or a index of nodes vector
}

bool BVH::Traverse(Ray& ray, Object** hit_obj, Vector& hit_point) {
	float tmp;
	float tmin = FLT_MAX; 
	bool hit = false;

	BVHNode* currentNode = nodes[0];

	*hit_obj = nullptr;

	// If the ray does not intersect the AABB of the root node return false
	if (!currentNode->getAABB().intercepts(ray, tmp)) {
		return false;
	}

	// Infinite loop to traverse the BVH tree until a leaf node is reached or the stack is empty
	for (;;) {
		if (!currentNode->isLeaf()) {

			float t_right = FLT_MAX, t_left = FLT_MAX;

			BVHNode* childLeft = nodes[currentNode->getIndex()];
			BVHNode* childRight = nodes[currentNode->getIndex() + 1];

			bool inter_left = childLeft->getAABB().intercepts(ray, t_left);
			bool inter_right = childRight->getAABB().intercepts(ray, t_right);


			StackItem* item = nullptr;

			// If both children intersect the ray, push the closest child becomes the current node 
			// and the farthest child is pushed to the stack
			if (inter_left && inter_right) {

				if (t_left > t_right) {
					currentNode = childRight;
					item = new StackItem(childLeft, t_left);
				}
				else {
					currentNode = childLeft;
					item = new StackItem(childRight, t_right);
				}
				hit_stack.push(*item);
				continue;
			}
			else if (inter_left) {
				currentNode = childLeft;
				continue;
			}
			else if (inter_right) {
				currentNode = childRight;
				continue;
			}
		}
		// If the current node is a leaf node, check for intersection with the objects in the leaf node
		else {
			int start = currentNode->getIndex(), objs = currentNode->getNObjs();
			for (int i = start; i < start + objs; i++) {
				if (objects.at(i)->intercepts(ray, tmp) && tmp < tmin) {
					*hit_obj = objects[i];
					tmin = tmp;
				}
			}
		}
	
		// Pop items from stack until it's empty
		bool end = false;
		for (;;) {
			if (hit_stack.empty()) {
				if (hit_obj == nullptr) return false;
				end = true;
				break;
			}
			else {
				StackItem item = hit_stack.top();
				hit_stack.pop();
				if (tmin > item.t) {
					currentNode = item.ptr;
					break;
				}
			}
		}
		if (end) break;
	}

	if (*hit_obj == nullptr) return false;

	hit_point = ray.origin + ray.direction * tmin;
	return true;
}

bool BVH::Traverse(Ray& ray) {  // Shadow ray with length
	float tmp;

	double length = ray.direction.length(); // Distance between light and intersection point
	ray.direction.normalize();

	BVHNode* currentNode = nodes[0];

	if (!currentNode->getAABB().intercepts(ray, tmp)) return false;

	// Infinite loop to traverse the BVH tree until a leaf node is reached or the stack is empty
	for (;;) {
		if (!currentNode->isLeaf()) {
			float t_right = FLT_MAX, t_left = FLT_MAX;

			BVHNode* childLeft = nodes[currentNode->getIndex()];
			BVHNode* childRight = nodes[currentNode->getIndex() + 1];

			bool inter_left = childLeft->getAABB().intercepts(ray, t_left);
			bool inter_right = childRight->getAABB().intercepts(ray, t_right);

			StackItem* item = nullptr;

			if (inter_left && inter_right) {

				currentNode = childLeft;
				item = new StackItem(childRight, t_right);
				hit_stack.push(*item);
				continue;

			}
			else if (inter_left) {

				currentNode = childLeft;
				continue;

			}
			else if (inter_right) {

				currentNode = childRight;
				continue;
			}
		}
		else {
			int start = currentNode->getIndex(), objs = currentNode->getNObjs();
			for (int i = start; i < start + objs; i++) {
				if (objects[i]->intercepts(ray, tmp) && tmp < length) {
					return true;
				}
			}
		}

		if (hit_stack.empty()) {
			return false;
		}
		else {
			StackItem item = hit_stack.top();
			hit_stack.pop();
			currentNode = item.ptr;
		}


	}
}
