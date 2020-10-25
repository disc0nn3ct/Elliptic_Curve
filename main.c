#include <gcrypt.h>
#include <stdio.h>
#include "e_curve.h"


int main()
{

	struct mong_curve m_c;

	set_parameters(&m_c);


	gcry_mpi_dump(m_c.p);
	printf(" p\n");
	gcry_mpi_dump(m_c.point1.Z);
	printf(" z\n");

	gcry_mpi_dump(m_c.point1.X);
	printf(" x\n");

	gcry_mpi_dump(m_c.point1.Y);
	printf(" y\n");

	gcry_mpi_dump(m_c.A);
	printf(" a\n");

	gcry_mpi_dump(m_c.B);
	printf(" a\n");



    printf("finished\n");
}