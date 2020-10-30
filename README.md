# Elliptic_Curve 

Вычислить кратную точку на эллиптической кривой в форме **Монтгомери**

Выполнено в Kali GNU/Linux Rolling 64 bit. 

Установка: [libgcrypt](https://gnupg.org/download/index.html) 

Нужно: Libgpg-error, Libgcrypt

Установка: сначала Libgpg-error -> Libgcrypt


 ```
./configurate 
make 
make check
make install 

 ```

Компиляция: gcc main.c e_curve.c -lgcrypt

main.c основной + тесты

e_curve все алгоритмы





