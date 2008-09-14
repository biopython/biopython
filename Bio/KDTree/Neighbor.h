/* Neighbor */

struct Neighbor;

struct Neighbor
{
        long int index1;
        long int index2;
        float radius;
	struct Neighbor* next;
};
