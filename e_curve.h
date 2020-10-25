#ifndef E_CURVE_H
#define E_CURVE_H


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

void init_point(struct point* p);

// инициализация переменных
void init_mong_curv(struct mong_curve* s);


// тут задается кривая и передаются параметры 

void set_parameters(struct mong_curve* m_c);

// Удвоение точки будет реализовано: http://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html
// Алгоритм dbl-1987-m-3 
void doubling_point(struct mong_curve* m_c);


// Функция для проверки точки на кривой.  Если точна находится на кривой, то возвращает 0, и если нет, то 1 

int is_point_on_curve(struct mong_curve* m_c);





#endif 