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

// создание копии
void make_copy_of_curve(struct mong_curve* copy, struct mong_curve* orig);


// Инициализация переменных для точки (для работы с большими числами) 

void init_point(struct point* p);

// инициализация переменных
void init_mong_curv(struct mong_curve* s);


// тут задается кривая и передаются параметры 

void set_parameters(struct mong_curve* m_c);


// Переход из проективных координат 
void transform_point(struct point* point_1, gcry_mpi_t* p);


// Удвоение точки будет реализовано: http://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html
// Алгоритм dbl-1987-m-3 
void doubling_point(struct mong_curve* m_c);


// Функция для проверки точки на кривой.  Если точна находится на кривой, то возвращает 0, и если нет, то 1 

int is_point_on_curve(struct mong_curve* m_c);


//dadd-1987-m  (Сложение), результат сложения в 1 переменной. p это модуль
void add_point(struct point* point_3, struct point* point_2, struct point* def, gcry_mpi_t* p );


// Реализация лестницы Монтгомери из лекций 53-56 слайды, реализация на 54 стр. 
// https://drive.google.com/drive/u/1/folders/17_F1NM91KR-6HOnUG_pHWxmtlRf7auU1	
void montgomery_ladder(struct point* point_1, gcry_mpi_t* k, struct mong_curve* m_c);



#endif 