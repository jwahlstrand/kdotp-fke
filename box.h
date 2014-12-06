#include <glib.h>
#include "three-vector.h"

/* represents a square box */

#ifndef BOX_H
#define BOX_H

enum {
    BOX_GOOD = 0,
    BOX_BAD_CENTER,
    BOX_BAD_EDGE
};

typedef struct {
    ThreeVector * center;
    ThreeVector * normal;
    double width;
} Box;

Box * new_box (const ThreeVector * center, const ThreeVector * normal, double width);

gboolean is_vector_inside_box (const ThreeVector * kp, const Box * b);
gboolean is_vector_on_edge_of_box (const ThreeVector * kp, const Box * b);
gboolean is_vector_on_corner_of_box (const ThreeVector * kp, const Box * box);
double box_area (const Box * box);

gint is_box_bad (const Box * b);

#endif

