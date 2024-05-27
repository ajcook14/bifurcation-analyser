
/*  error codes  */
#define NEWTON_DOUBLE_COUNT -1
#define NEWTON_EPSILON_WIDTH -2
#define ERROR_NO_BIFURCATION -3
#define ERROR_MAX_DERIVATIVE -4
#define ERROR_MAX_SUBDIVISIONS -5   // thrown when refine_measure reaches REFINE_MAX_SUBDIVISIONS
#define ERROR_AUTO_TOLERANCE -6     // thrown when automatic reaches AUTO_TOLERANCE

/*  parameters  */
#define REFINE_MAX_SUBDIVISIONS 4   // timeout limit when refine_measure isn't discarding any intervals
#define AUTO_TOLERANCE 1e-5         // minimum tolerance automatic is allowed to reach
