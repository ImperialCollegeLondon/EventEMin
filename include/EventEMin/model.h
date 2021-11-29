#ifndef EVENT_EMIN_MODEL_ALL_H
#define EVENT_EMIN_MODEL_ALL_H

// in 2D (image) space
#include "EventEMin/model/model/affinity.h"
#include "EventEMin/model/model/homography.h"
#include "EventEMin/model/model/isometry.h"
#include "EventEMin/model/model/rotation.h"
#include "EventEMin/model/model/similarity.h"
#include "EventEMin/model/model/translation.h"
#include "EventEMin/model/model/translation2d.h"
#include "EventEMin/model/model/translation_normal.h"
// in 2D (image) space incremental
#include "EventEMin/model/incremental_model/affinity.h"
#include "EventEMin/model/incremental_model/isometry.h"
#include "EventEMin/model/incremental_model/rotation.h"
#include "EventEMin/model/incremental_model/similarity.h"
#include "EventEMin/model/incremental_model/translation2d.h"
#include "EventEMin/model/incremental_model/translation_normal.h"

// in 3D space
#include "EventEMin/model/model/six_dof.h"
#include "EventEMin/model/model/translation3d.h"
// in 3D space incremental
#include "EventEMin/model/incremental_model/six_dof.h"
#include "EventEMin/model/incremental_model/translation3d.h"

#endif  // EVENT_EMIN_MODEL_ALL_H
