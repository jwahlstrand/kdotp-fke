#include <math.h>
#include <glib.h>
#include "box.h"

Box * new_box (const ThreeVector * k, const ThreeVector * normal, double width) {
	Box * b = g_slice_new (Box);
	b->center = three_vector_copy(k);
	b->normal = three_vector_copy(normal);
	g_assert(width>0);
	b->width = width;
	return b;
}

double box_area (const Box * box) {
    return (box->width)*(box->width);
}

gboolean is_vector_inside_box (const ThreeVector * kp, const Box * box) {
	gboolean by, bz;
	Box b = *box;
	by = ((kp->x[0]>(b.center->x[0]-b.width/2)) && (kp->x[0]<(b.center->x[0]+b.width/2)));
	bz = ((kp->x[1]>(b.center->x[1]-b.width/2)) && (kp->x[1]<(b.center->x[1]+b.width/2)));
	return (by && bz);
}

gboolean is_vector_on_edge_of_box (const ThreeVector * kp, const Box * box) {
    gboolean by, bz;
    Box b = *box;
    double edgep = b.center->x[0]+b.width/2;
    double edgem = b.center->x[0]-b.width/2;
    by = (fabs(kp->x[0]-edgep)<1e-6) || (fabs(kp->x[0]-edgem)<1e-6);
    edgep = b.center->x[1]+b.width/2;
    edgem = b.center->x[1]-b.width/2;
    bz = (fabs(kp->x[1]-edgep)<1e-6) || (fabs(kp->x[1]-edgem)<1e-6);
    return (by && bz);
}

gboolean is_vector_on_corner_of_box (const ThreeVector * kp, const Box * box) {
    double cx, cy;
    double kx = kp->x[0];
    double ky = kp->x[1];
    double w = box->width;
    cx=kx+w/2;
    cy=ky+w/2;
    if (fabs(kp->x[0]-cx)<1e-6 && fabs(kp->x[1]-cy)<1e-6)
        return TRUE;
    cx=kx+w/2;
    cy=ky-w/2;
    if (fabs(kp->x[0]-cx)<1e-6 && fabs(kp->x[1]-cy)<1e-6)
        return TRUE;
    cx=kx-w/2;
    cy=ky+w/2;
    if (fabs(kp->x[0]-cx)<1e-6 && fabs(kp->x[1]-cy)<1e-6)
        return TRUE;
    cx=kx-w/2;
    cy=ky-w/2;
    if (fabs(kp->x[0]-cx)<1e-6 && fabs(kp->x[1]-cy)<1e-6)
        return TRUE;
    return FALSE;
}

gint is_box_bad (const Box * b) {
    double kx = b->center->x[0];
    double ky = b->center->x[1];
    double w = b->width;
    if (fabs(fabs(kx)-fabs(ky))<1e-6)
        return BOX_BAD_CENTER;
    /*if ((fabs(kx-w/2)<1e-6) || (fabs(ky-w/2)<1e-6) || (fabs(kx+w/2)<1e-6) || (fabs(ky+w/2)<1e-6)) {
        return BOX_BAD_EDGE;
    }*/
    double cx, cy;
    cx=kx+w/2;
    cy=ky+w/2;
    if (fabs(fabs(cx)-fabs(cy))<1e-6)
        return BOX_BAD_EDGE;
    cx=kx+w/2;
    cy=ky-w/2;
    if (fabs(fabs(cx)-fabs(cy))<1e-6)
        return BOX_BAD_EDGE;
    cx=kx-w/2;
    cy=ky+w/2;
    if (fabs(fabs(cx)-fabs(cy))<1e-6)
        return BOX_BAD_EDGE;
    cx=kx-w/2;
    cy=ky-w/2;
    if (fabs(fabs(cx)-fabs(cy))<1e-6)
        return BOX_BAD_EDGE;
    return BOX_GOOD;
}
