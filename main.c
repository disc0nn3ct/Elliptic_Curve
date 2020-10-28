#include <gcrypt.h>
#include <stdio.h>
#include "e_curve.h"


int main()
{

	int flag_for_size = set_flag(0);  // Реализовано 2 режима: 0 это 256, а 1 это  512 

	struct mong_curve m_c;

	set_parameters(&m_c); // Установить гостовые переведенные параметры

	// Проверка переведенных данных из госта  (Начальные данные)
	printf("Тест 1. \n");
	printf("Проверка начальных значений\n");

	if(is_point_on_curve(&m_c) == 1 )
	{
		printf("Test Passed\n");
	}
	else
	{
		printf("Test Failed\n");
	}

	printf("\nТест 2.\n");
	printf("Проверяем на то, что [q]P = O, где q – порядок группы точек.\n");
	// Проверить, что [q]P = O, где q – порядок группы точек.
	set_parameters(&m_c); // Установить гостовые переведенные параметры



 	gcry_mpi_t q = gcry_mpi_new(0);
	if (flag_for_size == 0)
	{
		gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67", 0, 0);
	}
	else
	{
		gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC98CDBA46506AB004C33A9FF5147502CC8EDA9E7A769A12694623CEF47F023ED", 0, 0);
	}

	montgomery_ladder(&m_c.point1, &q, &m_c);

	gcry_mpi_t zero = gcry_mpi_new(0);
 	gcry_mpi_t one = gcry_mpi_new(0);


 	gcry_mpi_scan(&zero, GCRYMPI_FMT_HEX, "0", 0, 0);
 	gcry_mpi_scan(&one, GCRYMPI_FMT_HEX, "1", 0, 0);

	if( gcry_mpi_cmp(m_c.point1.Z, zero) ==0 && gcry_mpi_cmp(m_c.point1.X, one)==0 )
	{
		printf("Test Passed\n");
	}
	else
	{
		printf("Test Failed\n");
	}




	// print_point(&m_c.point1);







	set_parameters(&m_c); // Установить гостовые переведенные параметры








	printf("\nTest 3\n");










//////////////////////////////////////////////////////////////////////////////////////////////

	printf("\n\nТест 4.\n");

	// struct point point_k1;
	struct point point_k2;
	struct point point_k3;
	struct point point_k4;
	struct point point_k5;

	init_point(&point_k2);
	init_point(&point_k3);
	init_point(&point_k4);
	init_point(&point_k5);

	make_copy_point(&point_k2, &m_c.point1);
	make_copy_point(&point_k3, &m_c.point1);
	make_copy_point(&point_k4, &m_c.point1);
	make_copy_point(&point_k5, &m_c.point1);



	gcry_mpi_t k1 = gcry_mpi_new(0);
 	gcry_mpi_scan(&k1, GCRYMPI_FMT_HEX, "4", 0, 0);
	gcry_mpi_t k2 = gcry_mpi_new(0);
 	gcry_mpi_scan(&k2, GCRYMPI_FMT_HEX, "3", 0, 0);
	gcry_mpi_t k12 = gcry_mpi_new(0);
	gcry_mpi_addm(k12, k1, k2, m_c.p_mod);
	gcry_mpi_t summ = gcry_mpi_new(0);
	gcry_mpi_t subm = gcry_mpi_new(0);



	if(gcry_mpi_cmp(k1, k2) < 0) // Нужно, чтобы k1 было больше k2 
	{
		gcry_mpi_swap(k1,k2);
	}


	montgomery_ladder(&point_k2, &k1, &m_c);      										// [k1]P
	montgomery_ladder(&point_k3, &k2, &m_c);	  										// [k2]P
	gcry_mpi_addm(summ, k2, k1, m_c.p_mod);		  										// summ = k1 + k2 
	montgomery_ladder(&point_k4, &summ, &m_c); 	  										// [k1+k2]P
	gcry_mpi_subm(subm, k1,k2, m_c.p_mod); 		  										// [k1-k2]P   // так как формула сложения требует знание о еще одной точке 
	montgomery_ladder(&point_k5, &subm, &m_c);	  
	add_point(&point_k2, &point_k3, &point_k5, &m_c.p_mod);  							// [k1]P + [k2]P
	transform_point(&point_k2, &m_c.p_mod);


	if( (gcry_mpi_cmp(point_k2.X, point_k4.X) == 0)  && (gcry_mpi_cmp(point_k2.Z, point_k4.Z) == 0) )	
	{
		printf("Test Passed\n");
	}
	else
	{
		printf("Test Failed\n");
	}

	// Освобождаем память 
	del_point(&point_k2);
	del_point(&point_k3);
	del_point(&point_k4);
	del_point(&point_k5);



//////////////////////////////////////////////////////////////////////////////////////////
	// init_point(&point_k1);
	// init_point(&point_k2);
	// init_point(&point_k12);

	





	// montgomery_ladder(&point_k1, &k1, &m_c);

	// 	printf("\n\ngggggggggggggggg\n");

	// montgomery_ladder(&point_k2, &k2, &m_c);
	// printf("\n\ngggggggggggggggg\n");

	// montgomery_ladder(&point_k12, &k12, &m_c);


	// print_point(&point_k1);
	// print_point(&point_k2);
	
	// add_point(&summ, &point_k1, &point_k2, &m_c.p_mod);

	// print_point(&point_k12);



	// printf("\nTest 3\n");
	// make_copy_point(&point_k1, &m_c.point1);
	// montgomery_ladder(&point_k1, &m_c.p_mod, &m_c);
	// print_point(&point_k1);



	
	// montgomery_ladder();




	// gcry_mpi_dump(m_c.p_mod);
	// printf(" p\n");
	// gcry_mpi_dump(m_c.point1.Z);
	// printf(" z\n");

	// gcry_mpi_dump(m_c.point1.X);
	// printf(" x\n");

	// gcry_mpi_dump(m_c.point1.Y);
	// printf(" y\n");

	// gcry_mpi_dump(m_c.A);
	// printf(" A\n");

	// gcry_mpi_dump(m_c.B);
	// printf(" B\n");

	// printf(" \n\n\n\n\n");
 	// doubling_point(&m_c);
	// add_point(&m_c.point1, &m_c.point1, &m_c.point1, &m_c.p_mod);



 // 	gcry_mpi_dump(m_c.point1.Z);
	// printf(" z\n");

	// gcry_mpi_dump(m_c.point1.X);
	// printf(" x\n");

	// gcry_mpi_dump(m_c.point1.Y);
	// printf(" y\n");




	// is_point_on_curve(&m_c);

	// printf("ALLL NEW INFO +++++++++++++++++\n");


 // 	gcry_mpi_t four = gcry_mpi_new(0);
 // 	gcry_mpi_scan(&four, GCRYMPI_FMT_HEX, "4", 0, 0);
	// montgomery_ladder(&m_c.point1, &four, &m_c);


    printf("\n\nfinished\n");
}