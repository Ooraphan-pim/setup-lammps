
#ifndef DLO_MANAGER_STAFF_PARAMS_HH_LOADED
#define DLO_MANAGER_STAFF_PARAMS_HH_LOADED


const int  DLO_CHAR_BUF_SIZE = 256;
const int  DLO_MANAGER_RANK = 0;

// tags for messages from worker to dispatcher
const int  DLO_READY_TAG = 1;
const int  DLO_OUTPUT_TAG = 2;
const int  DLO_ERROR_ABORT_TAG = 3;

// tags for messages from dispatcher to worker
const int  DLO_ASSIGNMENT_TAG = 101;
const int  DLO_EXIT_TAG = 102;





#endif   // do not put anything after this #endif
