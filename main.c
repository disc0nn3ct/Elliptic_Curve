#include <gcrypt.h>
#include <stdio.h>
#include "e_curve.h"


int main()
{

	struct mong_curve m_c;

	set_parameters(&m_c);


	gcry_mpi_dump(m_c.p_mod);
	printf(" p\n");
	gcry_mpi_dump(m_c.point1.Z);
	printf(" z\n");

	gcry_mpi_dump(m_c.point1.X);
	printf(" x\n");

	gcry_mpi_dump(m_c.point1.Y);
	printf(" y\n");

	gcry_mpi_dump(m_c.A);
	printf(" A\n");

	gcry_mpi_dump(m_c.B);
	printf(" B\n");

	printf(" \n\n\n\n\n");
 	// doubling_point(&m_c);

 	gcry_mpi_dump(m_c.point1.Z);
	printf(" z\n");

	gcry_mpi_dump(m_c.point1.X);
	printf(" x\n");

	gcry_mpi_dump(m_c.point1.Y);
	printf(" y\n");




	is_point_on_curve(&m_c);

    printf("finished\n");
}