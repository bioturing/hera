#include "../attribute.h"

#define get_trans_len(t, i) ((t)->rlen[(i)])
#define get_trans_id(t, i) ((t)->rname + (t)->lname * (i))