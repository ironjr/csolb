#include "../core/solb.h"
#include "omp.h"

#define STRESS_TEST_NUM 100000000

void testSimple0();
void testStress0(int num);

int main()
{
	testStress0(STRESS_TEST_NUM);
	system("pause");
}

void testSimple0()
{
	top_solenoid_t *sol = new top_solenoid_t(.5, .5142, -.21, .21, 4.7619e8);
	//sol->print();
	mag_field_2d_t *field = solBSingle(sol, 0, 0);
	//field->print();
	delete field;
	delete sol;
}

void testStress0(int num)
{
	time_t begin;
	time_t end;
	struct tm *timeinfo;

	printf("Running total %d SolB queries.\n", num);

	time(&begin);
	timeinfo = localtime(&begin);
	printf("Test began at %s\n", asctime(timeinfo));

#pragma omp parallel for
    for (int i = 0; i < num; ++i)
    {
        testSimple0();
    }

	time(&end);
	timeinfo = localtime(&end);
	printf("Test ended at %s\n", asctime(timeinfo));

	double elapsed = difftime(end, begin);
	printf("\nSummary\n");
	printf("Time Elapsed         : %lf\n", elapsed);
	printf("Number of Queries    : %d\n", num);
	printf("Average Process Time : %lf\n", elapsed / num);
}
