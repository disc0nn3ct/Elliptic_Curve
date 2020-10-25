#include "e_curve.h"
#include <gcrypt.h>

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


	// Набор параметров id-tc26-gost-3410-2012-256-paramSetA (Будет реализовон переход из скрученной эллиптической кривой в форме Эдвардса в форму Монтгомери)
	gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97", 0, 0);
	gcry_mpi_scan(&e, GCRYMPI_FMT_HEX, "1", 0, 0);
	gcry_mpi_scan(&d, GCRYMPI_FMT_HEX, "605F6B7C183FA81578BC39CFAD518132B9DF62897009AF7E522C32D6DC7BFFB", 0, 0);
	gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67", 0, 0);
	gcry_mpi_scan(&u, GCRYMPI_FMT_HEX, "D", 0, 0);
	gcry_mpi_scan(&v, GCRYMPI_FMT_HEX, "60CA1E32AA475B348488C38FAB07649CE7EF8DBE87F22E81F92B2592DBA300E7", 0, 0);

	// Набор параметров id-tc26-g o s t-3410-2012-512-paramSetC  
	// gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7", 0, 0);
	// gcry_mpi_scan(&e, GCRYMPI_FMT_HEX, "1", 0, 0);
	// gcry_mpi_scan(&d, GCRYMPI_FMT_HEX, "9E4F5D8C017D8D9F13A5CF3CDF5BFE4DAB402D54198E31EBDE28A0621050439CA6B39E0A515C06B304E2CE43E79E369E91A0CFC2BC2A22B4CA302DBB33EE7550", 0, 0);
	// gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, "3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC98CDBA46506AB004C33A9FF5147502CC8EDA9E7A769A12694623CEF47F023ED", 0, 0);
	// gcry_mpi_scan(&u, GCRYMPI_FMT_HEX, "12", 0, 0);
	// gcry_mpi_scan(&v, GCRYMPI_FMT_HEX, "469AF79D1FB1F5E16B99592B77A01E2A0FDFB0D01794368D9A56117F7B38669522DD4B650CF789EEBF068C5D139732F0905622C04B2BAAE7600303EE73001A3D", 0, 0);

	// тут производилась отладка через gcry_mpi_dump() и проверка через wolfram mathematica. 

	m_c->p_mod=p;
	m_c->point1.Z = z; // z = 1

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
	// gcry_mpi_release(p); // так как используются, нельзя удалять 
	// gcry_mpi_release(z);				
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
 	

	printf("x=========================================\n");

	gcry_mpi_dump(m_c->p_mod);
	printf(" \nx=========================================\n");


 	gcry_mpi_addm(a_1, m_c->point1.X, m_c->point1.Z, m_c->p_mod); 	// X1 + Z1
  	gcry_mpi_mulm(a_1, a_1, a_1, m_c->p_mod);   					// (X1 + Z1)^2
	gcry_mpi_subm(b_1, m_c->point1.X, m_c->point1.Z, m_c->p_mod); 	// X1 - Z1
  	gcry_mpi_mulm(b_1, b_1, b_1, m_c->p_mod);   					// (X1 - Z1)^2
  	gcry_mpi_mulm(m_c->point1.X, a_1, b_1, m_c->p_mod);  			// X3 = (X1 + Z1)^2 * (X1 - Z1)^2
	gcry_mpi_subm(a_1, a_1, b_1 , m_c->p_mod);  		 			// (X1 + Z1)^2 - (X1 - Z1)^2
	gcry_mpi_invm(four, four,  m_c->p_mod); 			 			// 4^(-1)
	gcry_mpi_addm(a_24, m_c->A, two, m_c->p_mod);		 			// a + 2
	gcry_mpi_mulm(a_24, a_24, four, m_c->p_mod); 		 			//(a+2)/4
	gcry_mpi_mulm(a_24, a_24, a_1, m_c->p_mod);  		 			// (a+2)/4 * ((X1 + Z1)^2 - (X1 - Z1)^2)
	gcry_mpi_addm(a_24, b_1, a_24, m_c->p_mod); 		 			// (X1 - Z1)^2 + (a+2)/4 * ((X1 + Z1)^2 - (X1 - Z1)^2)
	gcry_mpi_mulm(m_c->point1.Z, a_1, a_24, m_c->p_mod); 			// ((X1 + Z1)^2 - (X1 - Z1)^2) * ((X1 - Z1)^2 + (a+2)/4 * ((X1 + Z1)^2 - (X1 - Z1)^2))

	gcry_mpi_release(a_1);
	gcry_mpi_release(b_1);	
	gcry_mpi_release(a_24);
	gcry_mpi_release(two);
	gcry_mpi_release(four);
}

// void add_point(struct mong_curve* m_c)
// {

// }



// Проверим точка находится на нашей кривой? Возвращает 0, если на кривой, 1 если нет. 
int is_point_on_curve(struct mong_curve* m_c)
{
	gcry_mpi_t l = gcry_mpi_new(0);
	gcry_mpi_t r = gcry_mpi_new(0);
	gcry_mpi_t r1 = gcry_mpi_new(0);

	gcry_mpi_mulm(l, m_c->point1.Y,m_c->point1.Y, m_c->p_mod); // y*y
	gcry_mpi_mulm(l, l, m_c->B, m_c->p_mod); // B*y^2
	gcry_mpi_mulm(r, m_c->point1.X, m_c->point1.X, m_c->p_mod); // x^2
	gcry_mpi_mulm(r, r, m_c->point1.X, m_c->p_mod); // x^3
	gcry_mpi_mulm(r1, m_c->point1.X, m_c->point1.X, m_c->p_mod); // x^2
	gcry_mpi_mulm(r1, r1, m_c->A, m_c->p_mod); // A*x^2
	gcry_mpi_addm(r, r, r1, m_c->p_mod); // x^3+A*x^2
	gcry_mpi_addm(r, r, m_c->point1.X, m_c->p_mod); // x^3+A*x^2+x
	
	if(gcry_mpi_cmp(l,r) == 0) 
	{
		printf("On our curve\n");
		return 0;
	}
	else
	{
		printf("Not in out curve\n");
		return 1;
	}
}










