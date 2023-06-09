# alignment_comparison
Программа для сравнения множественных выравниваний, созданная Артемом Тюкаевым, Яковом Коробицыным, Матвеем Киселевым и Георгием Малаховым. Алгоритм реализует сравнение двух множественных выравниваний одних и тех же последовательностей разными алгоритмами.

## Использование программы 
Программа предназначена для запуска из командной строки и принимает три обязательных аргумента:

`[путь до папки с программой]alignment_comparison.py [MA1] [MA2] [OUT]` - для командной строки Windows,

`python [путь до папки с программой]alignment_comparison.py [MA1] [MA2] [OUT]` - для командной строки Linux,

где [MA1] и [MA2] - пути до первого и второго выравниваний в формате FASTA соответственно, [OUT] - путь до файла вывода.

Файл вывода представляет собой tsv файл вида
```
F1-L1 F2-L2 - (шапка таблицы)
X1-Y1 X2-Y2
K1-M1 K2-M2
..... .....
```
где X1, K1 > 0 - координаты начал двух соседних блоков в первом выравнивании,
X2, K2 > 0 - координаты начал двух соседних блоков во втором выравнивании,
Y1, L1 > 0 - координаты концов двух соседних блоков в первом выравнивании,
Y2, L2 > 0 - координаты концов двух соседних блоков во втором выравнивании.

При запуске программы из командной строки Windows будьте внимательны: аргументы не должны быть в кавычках, а путь к Python должен быть одной из переменных среды.

Инструкции по поводу опций программы могут быть получены при вызове с опцией `-h`:

`alignment_comparison.py -h` - для командной строки Windows,

`python alignment_comparison.py -h` - для командной строки Linux.
