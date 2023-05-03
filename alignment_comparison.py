import bisect
import argparse
import textwrap as _textwrap


class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text_out = []
        for line in text.splitlines():
            for split_line in _textwrap.wrap(line, width):
                text_out.append(split_line)
        return text_out


parser = argparse.ArgumentParser(prog='alignment_comparison.py',
                                 usage='python %(prog)s [MA1] [MA2] [OUT]',
                                 formatter_class=LineWrapRawTextHelpFormatter,
                                 description='Алгоритм реализует сравнение двух множественных выравниваний одних и тех же последовательностей разными алгоритмами.\n'
                                             'Вывод представляет собой tsv файл из строк вида:\n'
                                             'F1-L1\tF2-L2 - (шапка таблицы)\n'
                                             'X1-Y1\tX2-Y2\n'
                                             'K1-M1\tK2-M2\n'
                                             '.....\t.....\n'
                                             'Где X1, K1 > 0 - координаты начал двух соседних блоков в первом выравнивании,\n'
                                             'X2, K2 > 0 - координаты начал двух соседних блоков во втором выравнивании,\n'
                                             'Y1, L1 > 0 - координаты концов двух соседних блоков в первом выравнивании,\n'
                                             'Y2, L2 > 0 - координаты концов двух соседних блоков во втором выравнивании.\n'
                                 )
parser.add_argument('MA1',
                    help='Путь до первого множественного выравнивания в формате FASTA')
parser.add_argument('MA2',
                    help='Путь до второго множественного выравнивания в формате FASTA')
parser.add_argument('OUT',
                    help='Путь до выходного файла')
parser.add_argument('-v',
                    '--visualise',
                    type=int,
                    nargs='?',
                    const=5,
                    help='Визуализация в человекочитаемом виде в STDOUT, демонстрирующая не отдельные совпадающие колонки, а целые блоки выравнивания.\n'
                         'VISUALISE - Количество блоков на одной строчке STDOUT, по умолчанию VISUALISE = 5\n'
                         'Выдача производится в формате:\n'
                         'X1-Y1\tK1-L1\t...\n'
                         'X2-Y2\tK2-L2\t...\n'
                         'Если какому-нибудь блоку в одном из выравниваний не соответствует блока в другом, то заместо координат блока во втором выравнивании ставятся прочерки.\n'
                    )
parser.add_argument('-p',
                    '--percent',
                    action='store_true',
                    help='Вывод в STDDOUT процент выравненных колонок в каждом из двух выравниваний.')

args = parser.parse_args()


def fasta_to_dist(fasta):
    seq_dict = {}
    name = ''
    for line in fasta:
        if line.startswith('>'):
            name = line[1: -1]
            seq_dict[name] = ''
        else:
            seq_dict[name] += line[: -1]

    for name, sequence in seq_dict.items():
        tokens = []
        for symbol in sequence:
            fol = -1 if len(tokens) == 0 else tokens[-1]
            if symbol != '-':
                tokens.append(fol + 1)
            else:
                tokens.append(fol)
        seq_dict[name] = tokens
    return seq_dict


def sort_dict(d):
    return sorted([i for i in d.keys()], reverse=True)


PATH1 = args.MA1
PATH2 = args.MA2

with open(PATH1) as alignment:
    first_alignment = fasta_to_dist(alignment)
with open(PATH2) as alignment:
    second_alignment = fasta_to_dist(alignment)

sorted_dict = sort_dict(first_alignment)
first_certainty = []
second_certainty = []


def f_symbol(seq, index):
    return first_alignment[sorted_dict[seq]][index]


def s_symbol(seq, index):
    return second_alignment[sorted_dict[seq]][index]


def compare_column(first_coordinate, second_coordinate):
    for seq in range(len(sorted_dict)):
        if f_symbol(seq, first_coordinate) == -1:
            f_gap = True
        elif first_coordinate != 0 and f_symbol(seq, first_coordinate) == f_symbol(seq, first_coordinate - 1):
            f_gap = True
        else:
            f_gap = False
        if s_symbol(seq, second_coordinate) == -1:
            s_gap = True
        elif second_coordinate != 0 and s_symbol(seq, second_coordinate) == s_symbol(seq, second_coordinate - 1):
            s_gap = True
        else:
            s_gap = False
        if f_symbol(seq, first_coordinate) != s_symbol(seq, second_coordinate) or f_gap ^ s_gap:
            first_certainty.append(-1)
            second_certainty.append(-1)
            return
    first_certainty.append(second_coordinate)
    second_certainty.append(first_coordinate)


def mark_alignment(seq, first_frame, second_frame):
    if len(first_frame) == 0 and len(second_frame) == 0:
        return
    if len(first_frame) == 0:
        for _ in range(len(second_frame)):
            second_certainty.append(-1)
        return
    if len(second_frame) == 0:
        for _ in range(len(first_frame)):
            first_certainty.append(-1)
        return

    if f_symbol(seq, first_frame[0]) == -1 or f_symbol(seq, first_frame[0]) == -1:
        subtraction = 0
    elif first_frame[0] != 0 and f_symbol(seq, first_frame[0]) == f_symbol(seq, first_frame[0] - 1):
        subtraction = f_symbol(seq, first_frame[0]) + 1
    elif second_frame[0] != 0 and s_symbol(seq, second_frame[0]) == s_symbol(seq, second_frame[0] - 1):
        subtraction = s_symbol(seq, second_frame[0]) + 1
    else:
        subtraction = 0

    if f_symbol(seq, first_frame[0]) > s_symbol(seq, second_frame[-1]):
        for _ in first_frame:
            first_certainty.append(-1)
        for _ in second_frame:
            second_certainty.append(-1)
        return
    if f_symbol(seq, first_frame[-1]) < s_symbol(seq, second_frame[0]):
        for _ in first_frame:
            first_certainty.append(-1)
        for _ in second_frame:
            second_certainty.append(-1)
        return

    if f_symbol(seq, first_frame[0]) > s_symbol(seq, second_frame[0]):
        shift = 0
        while s_symbol(seq, second_frame[shift]) != f_symbol(seq, first_frame[0]):
            shift += 1
        for _ in range(shift):
            second_certainty.append(-1)
        second_frame = second_frame[shift:]
    elif f_symbol(seq, first_frame[0]) < s_symbol(seq, second_frame[0]):
        shift = 0
        while s_symbol(seq, second_frame[0]) != f_symbol(seq, first_frame[shift]):
            shift += 1
        for _ in range(shift):
            first_certainty.append(-1)
        first_frame = first_frame[shift:]

        if f_symbol(seq, first_frame[-1]) == -1 + subtraction:
        if s_symbol(seq, second_frame[-1]) == -1 + subtraction:
            mark_alignment(seq + 1, first_frame, second_frame)
            return
        else:
            last_gap = bisect.bisect_left(second_alignment[sorted_dict[seq]], -1 + subtraction)
            mark_alignment(seq + 1, first_frame, second_frame[:last_gap])
            for _ in range(len(second_frame[last_gap:])):
                second_certainty.append(-1)
            return
    if s_symbol(seq, second_frame[-1]) == -1 + subtraction:
        last_gap = bisect.bisect_left(first_alignment[sorted_dict[seq]], -1 + subtraction)
        mark_alignment(seq + 1, first_frame[:last_gap], second_frame)
        for _ in range(len(first_frame[last_gap:])):
            first_certainty.append(-1)
        return

    if f_symbol(seq, 0) == -1 + subtraction:
        first_shift = 0
        while f_symbol(seq, first_frame[first_shift]) == -1 + subtraction:
            first_shift += 1
        second_shift = 0
        while s_symbol(seq, second_frame[second_shift]) == -1 + subtraction:
            second_shift += 1
        mark_alignment(seq + 1, first_frame[:first_shift], second_frame[:second_shift])
        first_frame = first_frame[first_shift:]
        second_frame = second_frame[second_shift:]

    compare_column(first_frame[0], second_frame[0])
    mark_alignment(seq, first_frame[1:], second_frame[1:])


mark_alignment(0, list(range(len(first_alignment[sorted_dict[0]]))), list(range(len(second_alignment[sorted_dict[0]]))))
output = []
for first, second in enumerate(first_certainty):
    if second != -1:
        output.append([first, second])


def count_clusters():
    f_clusters = []
    for i, match in enumerate(first_certainty):
        if match != -1:
            if len(f_clusters) == 0 or f_clusters[-1][2]:
                f_clusters.append([i + 1, i + 1, False])
            else:
                f_clusters[-1][1] += 1
        else:
            if len(f_clusters) == 0 or not f_clusters[-1][2]:
                f_clusters.append([i + 1, i + 1, True])
            else:
                f_clusters[-1][1] += 1

    s_clusters = []
    for i, match in enumerate(second_certainty):
        if match != -1:
            if len(s_clusters) == 0 or s_clusters[-1][2]:
                s_clusters.append([i + 1, i + 1, False])
            else:
                s_clusters[-1][1] += 1
        else:
            if len(s_clusters) == 0 or not s_clusters[-1][2]:
                s_clusters.append([i + 1, i + 1, True])
            else:
                s_clusters[-1][1] += 1
    return f_clusters, s_clusters


def visualise():
    f_clusters, s_clusters = count_clusters()
    f_print = []
    s_print = []
    for i in range(len(f_clusters)):
        if f_clusters[i][2]:
            out = f'{f_clusters[i][0]}-{f_clusters[i][1]}'
            f_print.append(out)
            s_print.append('-' * len(out))
            out = f'{s_clusters[i][0]}-{s_clusters[i][1]}'
            f_print.append('-' * len(out))
            s_print.append(out)
        else:
            f_out = f'{f_clusters[i][0]}-{f_clusters[i][1]}'
            s_out = f'{s_clusters[i][0]}-{s_clusters[i][1]}'
            if len(f_out) > len(s_out):
                s_out = s_out + ' ' * (len(f_out) - len(s_out))
            elif len(f_out) < len(s_out):
                f_out = f_out + ' ' * (len(s_out) - len(f_out))
            f_print.append(f_out)
            s_print.append(s_out)

    for blocks in range(0, len(f_print), args.visualise):
        print(*f_print[blocks: blocks + args.visualise], sep='\t')
        print(*s_print[blocks: blocks + args.visualise], sep='\t')
        print('')


def save(file):
    file.write('F1-L1\tF2-L2\n')
    f_clusters, s_clusters = count_clusters()
    for i in range(len(f_clusters)):
        if not f_clusters[i][2]:
            file.write(f'{f_clusters[i][0]}-{f_clusters[i][1]}\t{s_clusters[i][0]}-{s_clusters[i][1]}\n')


if args.visualise is not None:
    visualise()
if args.percent:
    print(f'Процент выравненных колонок в первом выравнивании - {len(output) / len(first_alignment[sorted_dict[0]]) * 100:.2f}')
    print(f'Процент выравненных колонок во втором выравнивании - {len(output) / len(second_alignment[sorted_dict[0]]) * 100:.2f}')

with open(args.OUT, 'w') as file:
    save(file)
