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

	set_parameters(&m_c); // Установить гостовые переведенные параметры, чтобы тесты были независимыми 

	printf("\nTest 3\n");
	printf("Проверить, что [q + 1]P = P и [q − 1]P = −P.\n");

	struct point proint_q1, proint_q2; 
	
	init_point(&proint_q1);
	init_point(&proint_q2);
	make_copy_point(&proint_q1, &m_c.point1);
	make_copy_point(&proint_q2, &m_c.point1);


 	gcry_mpi_t q1 = gcry_mpi_new(0);
 	gcry_mpi_t q2 = gcry_mpi_new(0);
 	q1 = gcry_mpi_copy(q);
 	q2 = gcry_mpi_copy(q);

 	gcry_mpi_addm(q1, q1, one, m_c.p_mod); 						// q + 1
 	gcry_mpi_subm(q2, q2, one, m_c.p_mod); 						// q - 1

 	montgomery_ladder(&proint_q1, &q1, &m_c);
	montgomery_ladder(&proint_q2, &q2, &m_c);


	if(gcry_mpi_cmp(proint_q1.X, m_c.point1.X) == 0 && gcry_mpi_cmp(proint_q1.Y, m_c.point1.Y) == 0)
	{
		printf("Test Passed for [q+1]P = P\n");
	}
	else
	{
		printf("Test Failed\n");
	}
	

	if(gcry_mpi_cmp(proint_q2.X, m_c.point1.X) == 0 && gcry_mpi_cmp(proint_q2.Y, m_c.point1.Y) == 0)
	{
		printf("Test Passed for [q-1]P = -P\n");
	}
	else
	{
		printf("Test Failed\n");
	}





 	// очистка памяти
 	gcry_mpi_release(zero);
 	gcry_mpi_release(one);
 	gcry_mpi_release(q);
	gcry_mpi_release(q1);
	gcry_mpi_release(q2);
 	set_parameters(&m_c); // Установить гостовые переведенные параметры, чтобы тесты были независимыми 



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
	gcry_mpi_randomize(k1, 35, GCRY_STRONG_RANDOM); // рандомим k1 

 	// gcry_mpi_scan(&k1, GCRYMPI_FMT_HEX, "4", 0, 0);
 	
	gcry_mpi_t k2 = gcry_mpi_new(0);
	gcry_mpi_randomize(k2, 35, GCRY_STRONG_RANDOM);  // рандомим k2 

 	// gcry_mpi_scan(&k2, GCRYMPI_FMT_HEX, "3", 0, 0);
	
	gcry_mpi_t k12 = gcry_mpi_new(0);
	gcry_mpi_addm(k12, k1, k2, m_c.p_mod);
	gcry_mpi_t summ = gcry_mpi_new(0);
	gcry_mpi_t subm = gcry_mpi_new(0);



	if(gcry_mpi_cmp(k1, k2) < 0) // Нужно, чтобы k1 было больше k2 
	{
		gcry_mpi_swap(k1,k2);
	}


	montgomery_ladder(&point_k2, &k1, &m_c);                    // [k1]P
	montgomery_ladder(&point_k3, &k2, &m_c);                    // [k2]P
	gcry_mpi_addm(summ, k2, k1, m_c.p_mod);                     // summ = k1 + k2 
	montgomery_ladder(&point_k4, &summ, &m_c);                  // [k1+k2]P
	gcry_mpi_subm(subm, k1,k2, m_c.p_mod);                      // [k1-k2]P   // так как формула сложения требует знание о еще одной точке 
	montgomery_ladder(&point_k5, &subm, &m_c);	  
	add_point(&point_k2, &point_k3, &point_k5, &m_c.p_mod);     // [k1]P + [k2]P
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

	printf("\nТест 5\n");
	printf("\nЛежит ли [k]P на кривой?\n");

 	gcry_mpi_t rand1 = gcry_mpi_new(0);
	gcry_mpi_randomize(rand1, 256, GCRY_STRONG_RANDOM);  // рандомим k2 

 	montgomery_ladder(&m_c.point1, &rand1, &m_c);



	is_point_on_curve(&m_c);

	// gcry_mpi_dump(rand1);


 	del_curve(&m_c);
	gcry_mpi_release(rand1);


	gcry_mpi_release(k1);
	gcry_mpi_release(k12);
	gcry_mpi_release(summ);
	gcry_mpi_release(subm);
	del_point(&proint_q1);
	del_point(&proint_q2);
	
	gcry_control(GCRYCTL_FINALIZE);

    printf("\n\nfinished\n");
}