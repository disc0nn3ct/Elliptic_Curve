#include "e_curve.h"
#include <gcrypt.h>


int flag_for_size = 0; // Переменная для выбора режима: 0 это 256, а 1 это 512 

int set_flag(int a)
{
	flag_for_size =a;
	return flag_for_size;
}

// структура для хранения точки
struct point
{
	gcry_mpi_t X; 
	gcry_mpi_t Y;  
	gcry_mpi_t Z;
};

// структура для хранения параметров кривой Монтгомери и точек порождающих ее 

struct mong_curve
{
	gcry_mpi_t A;
	gcry_mpi_t B;
	gcry_mpi_t p_mod;
	struct point point1; 
};

// Инициализация переменных для точки (для работы с большими числами) 

void init_point(struct point* p1)
{
	p1->X=gcry_mpi_new(0);
	p1->Y=gcry_mpi_new(0);
	p1->Z=gcry_mpi_new(0);
}
// инициализация переменных
void init_mong_curv(struct mong_curve* s)
{
	s->A = gcry_mpi_new(0);
	s->B = gcry_mpi_new(0);
	s->p_mod=gcry_mpi_new(0);		
	init_point(&s->point1);
}

// создание копии точки 
void make_copy_point(struct point* where, struct point* from)
{
	where->X = gcry_mpi_copy(from->X);
	where->Y = gcry_mpi_copy(from->Y);
	where->Z = gcry_mpi_copy(from->Z);
}

void make_copy_of_curve(struct mong_curve* copy, struct mong_curve* orig)
{
	copy->A = gcry_mpi_copy(orig->A);
	copy->B = gcry_mpi_copy(orig->B);
	copy->p_mod= gcry_mpi_copy(orig->p_mod); 
	copy->point1.X =  gcry_mpi_copy(orig->point1.X); 
	copy->point1.Y = gcry_mpi_copy(orig->point1.Y); 
	copy->point1.Z = gcry_mpi_copy(orig->point1.Z);
}

// Вывод всех координат точки 
void print_point(struct point* point_1)
{	
	printf(" \nx+++++++++++++++++++++++++++++++++\n");
	gcry_mpi_dump(point_1->X);
	printf(" X\n");
	gcry_mpi_dump(point_1->Y); // Так как отбрасываем, сделано для общего случая 
	printf(" Y\n");
	gcry_mpi_dump(point_1->Z);
	printf(" Z\n");	
	printf(" \nx+++++++++++++++++++++++++++++++++\n");
}


// Освобождение памяти выделенной для point
void del_point(struct point* point_1)
{
	gcry_mpi_release(point_1->X);
	gcry_mpi_release(point_1->Y);
	gcry_mpi_release(point_1->Z);
}

// Освобождение памяти выделенной для структуры кривой + точки
void del_curve(struct mong_curve* m_c)
{
	gcry_mpi_release(m_c->A);
	gcry_mpi_release(m_c->B);
	gcry_mpi_release(m_c->p_mod);
	del_point(&m_c->point1);
}



// Набор параметров id-tc26-gost-3410-2012-512-paramSet
//https://docplayer.ru/46408167-Tehnicheskiy-komitet-026-zadanie-parametrov-skruchennyh-ellipticheskih-krivyh-edvardsa-v-sootvetstvii-s-gost-r.html 
// зададим параметры 

void set_parameters(struct mong_curve* m_c)
{
	gcry_mpi_t buf = gcry_mpi_new(0);
	gcry_mpi_t p = gcry_mpi_new(0);
 	gcry_mpi_t e = gcry_mpi_new(0);
 	gcry_mpi_t d = gcry_mpi_new(0);
 	gcry_mpi_t q = gcry_mpi_new(0);
 	gcry_mpi_t u = gcry_mpi_new(0);
 	gcry_mpi_t v = gcry_mpi_new(0);
 	gcry_mpi_t z = gcry_mpi_new(0);
 	gcry_mpi_t one = gcry_mpi_new(0);
 	gcry_mpi_t two = gcry_mpi_new(0);
 	gcry_mpi_t four = gcry_mpi_new(0);
 	init_mong_curv(m_c);


 	gcry_mpi_scan(&one, GCRYMPI_FMT_HEX, "1", 0, 0);
 	gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, 0);
 	gcry_mpi_scan(&four, GCRYMPI_FMT_HEX, "4", 0, 0);
	gcry_mpi_scan(&z, GCRYMPI_FMT_HEX, "1", 0, 0);

	

	if (flag_for_size == 0)
	{
	// Набор параметров id-tc26-gost-3410-2012-256-paramSetA (Будет реализовон переход из скрученной эллиптической кривой в форме Эдвардса в форму Монтгомери)
		gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97", 0, 0);
		gcry_mpi_scan(&e, GCRYMPI_FMT_HEX, "1", 0, 0);
		gcry_mpi_scan(&d, GCRYMPI_FMT_HEX, "605F6B7C183FA81578BC39CFAD518132B9DF62897009AF7E522C32D6DC7BFFB", 0, 0);
		gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67", 0, 0);
		gcry_mpi_scan(&u, GCRYMPI_FMT_HEX, "D", 0, 0);
		gcry_mpi_scan(&v, GCRYMPI_FMT_HEX, "60CA1E32AA475B348488C38FAB07649CE7EF8DBE87F22E81F92B2592DBA300E7", 0, 0);
	}
	else
	{

	// Набор параметров id-tc26-g o s t-3410-2012-512-paramSetC  
		gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7", 0, 0);
		gcry_mpi_scan(&e, GCRYMPI_FMT_HEX, "1", 0, 0);
		gcry_mpi_scan(&d, GCRYMPI_FMT_HEX, "9E4F5D8C017D8D9F13A5CF3CDF5BFE4DAB402D54198E31EBDE28A0621050439CA6B39E0A515C06B304E2CE43E79E369E91A0CFC2BC2A22B4CA302DBB33EE7550", 0, 0);
		gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC98CDBA46506AB004C33A9FF5147502CC8EDA9E7A769A12694623CEF47F023ED", 0, 0);
		gcry_mpi_scan(&u, GCRYMPI_FMT_HEX, "12", 0, 0);
		gcry_mpi_scan(&v, GCRYMPI_FMT_HEX, "469AF79D1FB1F5E16B99592B77A01E2A0FDFB0D01794368D9A56117F7B38669522DD4B650CF789EEBF068C5D139732F0905622C04B2BAAE7600303EE73001A3D", 0, 0);
	}
	// тут производилась отладка через gcry_mpi_dump() и проверка через wolfram mathematica. 

	m_c->p_mod=gcry_mpi_copy(p);
	m_c->point1.Z = gcry_mpi_copy(z); // z = 1

	// Реализация перехода: 
	gcry_mpi_addm(m_c->point1.X, one, v, p);   								// 1 + v 
	gcry_mpi_subm(buf, one, v, p);   										// 1 - v 
	gcry_mpi_invm(buf, buf, p); 											// (1-v)^(-1)
	gcry_mpi_mulm(m_c->point1.X, m_c->point1.X, buf, p); 					// (1+v)/(1-v)

	gcry_mpi_addm(m_c->point1.Y, one, v, p);   								// 1 + v 
	gcry_mpi_subm(buf, one, v, p);   										// 1 - v 
	gcry_mpi_mulm(buf, buf, u, p);	 										// u*(1 -v)
	gcry_mpi_invm(buf, buf, p);  											// (u*(1 -v))^(-1)
	gcry_mpi_mulm(m_c->point1.Y, m_c->point1.Y, buf, p);

	gcry_mpi_addm(m_c->A, e, d, p);  										// e + d 
 	gcry_mpi_mulm(m_c->A, two, m_c->A, p); 									// 2*(e + d)
	gcry_mpi_subm(buf, e, d, p); 											// e - d 
	gcry_mpi_invm(buf, buf, p); 											// (e-d)^(-1)
	gcry_mpi_mulm(m_c->A, m_c->A, buf, p);								    // 2*(e + d)/(e-d) 
 
	gcry_mpi_mulm(m_c->B, four, buf, p); 									// 4/(e-d)

	// Удалеие(освобождение) выделеной памяти

	gcry_mpi_release(buf);
	gcry_mpi_release(e);
	gcry_mpi_release(d);
	gcry_mpi_release(q);
	gcry_mpi_release(u);
	gcry_mpi_release(v);
	gcry_mpi_release(one);
	gcry_mpi_release(two);
	gcry_mpi_release(four);
	gcry_mpi_release(p);  
	gcry_mpi_release(z);				
}

// Переход из проективных координат 
void transform_point(struct point* point_1, gcry_mpi_t* p)
{
	gcry_mpi_t littl = gcry_mpi_new(0);
	gcry_mpi_t zeroo = gcry_mpi_new(0);
	gcry_mpi_scan(&zeroo, GCRYMPI_FMT_HEX, "0", 0, 0);

	if(gcry_mpi_cmp(point_1->Z, zeroo) != 0)
	{	
		// gcry_mpi_abs( point_1->Z);
		gcry_mpi_invm(littl, point_1->Z, *p);
		gcry_mpi_mulm(point_1->X, littl, point_1->X, *p);
		gcry_mpi_mulm(point_1->Z, littl, point_1->Z, *p);
	}
	else
	{
		gcry_mpi_invm(littl, point_1->X, *p);
		gcry_mpi_mulm(point_1->X, littl, point_1->X, *p);
	}

}


// Удвоение точки будет реализовано: http://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html
// Алгоритм dbl-1987-m-3
void doubling_point(struct mong_curve* m_c)
{
 	gcry_mpi_t a_1 = gcry_mpi_new(0);
 	gcry_mpi_t b_1 = gcry_mpi_new(0);
	gcry_mpi_t two = gcry_mpi_new(0);
 	gcry_mpi_t four = gcry_mpi_new(0);
 	gcry_mpi_t a_24 = gcry_mpi_new(0);

 	gcry_mpi_scan(&two, GCRYMPI_FMT_HEX, "2", 0, 0);
 	gcry_mpi_scan(&four, GCRYMPI_FMT_HEX, "4", 0, 0);
 	


 	gcry_mpi_addm(a_1, m_c->point1.X, m_c->point1.Z, m_c->p_mod); 	// X1 + Z1
  	gcry_mpi_mulm(a_1, a_1, a_1, m_c->p_mod);   					// (X1 + Z1)^2
	gcry_mpi_subm(b_1, m_c->point1.X, m_c->point1.Z, m_c->p_mod); 	// X1 - Z1
  	gcry_mpi_mulm(b_1, b_1, b_1, m_c->p_mod);   					// (X1 - Z1)^2
  	gcry_mpi_mulm(m_c->point1.X, a_1, b_1, m_c->p_mod);  			// X3 = (X1 + Z1)^2 * (X1 - Z1)^2
	gcry_mpi_subm(a_1, a_1, b_1, m_c->p_mod);  		 				// (X1 + Z1)^2 - (X1 - Z1)^2

	gcry_mpi_invm(four, four,  m_c->p_mod); 			 			// 4^(-1) 
	gcry_mpi_add(a_24, m_c->A, two);		 						// a + 2
	gcry_mpi_mulm(a_24, a_24, four, m_c->p_mod); 					// (a+2)/4    //  Модуль обязательно! ( проверено через wolfram mathematica )

	gcry_mpi_mulm(a_24, a_24, a_1, m_c->p_mod);  		 			// (a+2)/4 * ((X1 + Z1)^2 - (X1 - Z1)^2)
	gcry_mpi_addm(a_24, b_1, a_24, m_c->p_mod); 		 			// (X1 - Z1)^2 + (a+2)/4 * ((X1 + Z1)^2 - (X1 - Z1)^2)
	gcry_mpi_mulm(m_c->point1.Z, a_1, a_24, m_c->p_mod); 			// ((X1 + Z1)^2 - (X1 - Z1)^2) * ((X1 - Z1)^2 + (a+2)/4 * ((X1 + Z1)^2 - (X1 - Z1)^2))

	gcry_mpi_release(a_1);
	gcry_mpi_release(b_1);	
	gcry_mpi_release(a_24);
	gcry_mpi_release(two);
	gcry_mpi_release(four);
}


// Проверим точка находится на нашей кривой? Возвращает 0, если на кривой, 1 если нет. 
int is_point_on_curve(struct mong_curve* m_c)
{
	gcry_mpi_t l = gcry_mpi_new(0);
	gcry_mpi_t r = gcry_mpi_new(0);
	gcry_mpi_t r1 = gcry_mpi_new(0);

	gcry_mpi_mulm(l, m_c->point1.Y,m_c->point1.Y, m_c->p_mod); 						// y*y
	gcry_mpi_mulm(l, l, m_c->B, m_c->p_mod); 										// B*y^2
	gcry_mpi_mulm(r, m_c->point1.X, m_c->point1.X, m_c->p_mod); 					// x^2
	gcry_mpi_mulm(r, r, m_c->point1.X, m_c->p_mod); 								// x^3
	gcry_mpi_mulm(r1, m_c->point1.X, m_c->point1.X, m_c->p_mod); 					// x^2
	gcry_mpi_mulm(r1, r1, m_c->A, m_c->p_mod); 										// A*x^2
	gcry_mpi_addm(r, r, r1, m_c->p_mod); 											// x^3+A*x^2
	gcry_mpi_addm(r, r, m_c->point1.X, m_c->p_mod);								    // x^3+A*x^2+x
	
	if(gcry_mpi_cmp(l,r) == 0) 
	{
		printf("On our curve\n");
		return 1;
	}
	else
	{
		printf("Not in out curve\n");
		return 0;
	}
}



//dadd-1987-m-2 (Сложение) // результат будет записан в point_3 
void add_point(struct point* point_3, struct point* point_2, struct point* def, gcry_mpi_t* p )
{
	// Для тестов реализую 2 алгоритма http://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html
	int a = 0;  // это флаг для выбора релизации суммы 

 	gcry_mpi_t a_1 = gcry_mpi_new(0);
 	gcry_mpi_t a_2 = gcry_mpi_new(0);
 	gcry_mpi_t a_3 = gcry_mpi_new(0);
 	gcry_mpi_t a_4 = gcry_mpi_new(0);

	if(a == 0)
	{
	//dadd-1987-m-2 

	 	gcry_mpi_subm(a_1, point_3->X, point_3->Z, *p); 									// X3 - Z3
	 	gcry_mpi_addm(a_2, point_2->X, point_2->Z, *p); 									// X2 + Z2
	 	gcry_mpi_addm(a_3, point_3->X, point_3->Z, *p); 									// X3 + Z3
	 	gcry_mpi_subm(a_4, point_2->X, point_2->Z, *p); 									// X2 - Z2
	 	gcry_mpi_mulm(a_1, a_1, a_2, *p); 													// (X3 - Z3)*(X2 + Z2)
		gcry_mpi_mulm(a_3, a_3, a_4, *p); 													// (X3 - Z3)*(X2 + Z2)

		gcry_mpi_addm(a_2, a_1, a_3, *p); 													// (X3 - Z3)*(X2 + Z2) + (X3 - Z3)*(X2 + Z2)
		gcry_mpi_mulm(a_2, a_2, a_2, *p);  													// ((X3 - Z3)*(X2 + Z2) + (X3 - Z3)*(X2 + Z2))^2
		gcry_mpi_mulm(point_3->X, def->Z, a_2, *p);											// Z1 * ((X3 - Z3)*(X2 + Z2) + (X3 - Z3)*(X2 + Z2))^2
		
		gcry_mpi_subm(a_2, a_1, a_3, *p); 													// (X3 - Z3)*(X2 + Z2) + (X3 - Z3)*(X2 + Z2)
		gcry_mpi_mulm(a_2, a_2, a_2, *p);  													// ((X3 - Z3)*(X2 + Z2) + (X3 - Z3)*(X2 + Z2))^2
		gcry_mpi_mulm(point_3->Z, def->X, a_2, *p); 										// Z1 * ((X3 - Z3)*(X2 + Z2) + (X3 - Z3)*(X2 + Z2))^2
		
 	} 
 	else
 	{
 		//dadd-1987-m
 		printf("XXXXXXXXXXXXXXXX22222222222222222\n");
 		print_point(point_2);
 		printf("XXXXXXXXXXXXXXXX33333333333333333\n");
 		print_point(point_3);
 		gcry_mpi_mulm(a_1, point_2->X, point_3->X, *p);  									// X2*X3
		
 		printf("fffffffffffffffffffffffff\n");
 		gcry_mpi_dump(a_1);
		printf("\nfffffffffffffffffffffffff\n");
		gcry_mpi_mulm(a_2, point_2->Z, point_3->Z, *p); 								    // Z2*Z3
		gcry_mpi_subm(a_1, a_1, a_2, *p); 													// X2*X3 - Z2*Z3
		gcry_mpi_mulm(a_1, a_1, a_1, *p); 													// (X2*X3 - Z2*Z3)^2
		gcry_mpi_mulm(a_4, def->Z, a_1, *p); 												// Z0 * (X2*X3 - Z2*Z3)^2

		gcry_mpi_mulm(a_1, point_2->X, point_3->Z, *p);  									// X2*Z3
		gcry_mpi_mulm(a_2, point_2->Z, point_3->X, *p);  									// Z2*X3
		gcry_mpi_subm(a_1, a_1, a_2, *p); 													// X2*Z3 - Z2*X3		
		gcry_mpi_mulm(a_1, a_1, a_1, *p); 													// (X2*Z3 - Z2*X3)^2
		gcry_mpi_mulm(point_3->Z, def->X , a_1, *p); 										//  X1 * (X2*Z3 - Z2*X3)^2
		point_3->X = gcry_mpi_copy(a_4);
 	}

  // освобождение памяти
  	gcry_mpi_release(a_1);
  	gcry_mpi_release(a_2);
  	gcry_mpi_release(a_3);
  	gcry_mpi_release(a_4);
}

// Реализация лестницы Монтгомери из лекций 53-56 слайды, реализация на 54 стр. 
// https://drive.google.com/drive/u/1/folders/17_F1NM91KR-6HOnUG_pHWxmtlRf7auU1
	
void montgomery_ladder(struct point* point_old, gcry_mpi_t* k, struct mong_curve* m_c)
{
	unsigned int bitrate = gcry_mpi_get_nbits(*k);  // получаем длину бинарного вида
	struct mong_curve point_new;
	init_mong_curv(&point_new);

	struct mong_curve point0;          
	init_mong_curv(&point0);

	make_copy_of_curve(&point_new, m_c);  
	make_copy_of_curve(&point0, m_c);


	gcry_mpi_scan(&point0.point1.X, GCRYMPI_FMT_HEX, "1", 0, 0);
	gcry_mpi_scan(&point0.point1.Z, GCRYMPI_FMT_HEX, "0", 0, 0);

	for(int i=bitrate-1; i>=0; i--)
	{
		if(gcry_mpi_test_bit(*k, i) == 1)   
		{
			add_point(&point0.point1, &point_new.point1, point_old, &m_c->p_mod);
			doubling_point(&point_new);
		}
		else
		{
			add_point(&point_new.point1, &point0.point1, point_old, &m_c->p_mod);
			doubling_point(&point0);
		}
	}

	point_old->X =gcry_mpi_copy(point0.point1.X);
	point_old->Z =gcry_mpi_copy(point0.point1.Z);

	transform_point(point_old, &m_c->p_mod);


	// очищаем память 
	del_curve(&point_new);	
	del_curve(&point0);	
}






