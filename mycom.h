#define EC_No_errors            0
#define EC_Small_param_count    1
#define EC_Bad_file_name        2
#define EC_Many_open_files      3
#define EC_Bad_file_access      4
#define EC_File_reading_error   5
#define EC_File_writing_error   6
#define EC_Close_file_error     7
#define EC_Fatal_error          8
#define EC_No_enough_memory     9
#define EC_Bad_data_format     10
#define EC_No_data             11

#define pi 3.14159265358979323846264338327950280

int iabs(int a);

int imax(int a, int b);

int imin(int a, int b);

double dabs(double a);

double dmax(double a, double b);

double dmin(double a, double b);

double dsin(double x);

double dcos(double x);

double dexp(double x);

double mytime(const int n);

void myerr(char *msg, const int n);

double integrate(double f(double x), double a, double b, int n);

int fopen_m(FILE** F, const char* name, const char* mode);

int fclose_m(FILE** F);
