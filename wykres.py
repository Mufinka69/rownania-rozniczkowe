import pandas as pd
import matplotlib.pyplot as plt


def open_file(plik, color = None, s = "s", alpha = 1, name = None, markersize = None):
    t = []
    y = []

    df = pd.read_csv(plik, sep=';', skipinitialspace=True)
    print(df.columns.tolist())
    if 't' not in df.columns:
        raise ValueError("Brakuje kolumny 't' (czas) w pliku.")
    for col in df.columns:
        if col != 't' and col.startswith('y'):
            if(col == "y0"):
                # plt.plot(df['t'], df[col], s , label = fr"$y_{{\text{{{col[1]}}}}}(t)$", alpha = alpha)
                # plt.plot(df['t'], df[col], s , label = fr"$y_{{\text{{{col[1]}}}}}$ - {name}", alpha = alpha, markersize = markersize)
                plt.plot(df['t'], df[col], s , label = fr"{name}", alpha = alpha, markersize = markersize, color = color)


def open_sol(plik, color = None, s = "s", alpha = 1, name = None, markersize = None):
    t = []
    y = []

    df = pd.read_csv(plik, sep=';', skipinitialspace=True)
    print(df.columns.tolist())
    if 't' not in df.columns:
        raise ValueError("Brakuje kolumny 't' (czas) w pliku.")
    for col in df.columns:
        if col != 't' and col.startswith('y'):
            # $y_{{\text{{{col[1]}}}}}$ - 
            # plt.plot(df['t'], df[col], s , label = fr"$y_{{\text{{{col[1]}}}}}(t)$", alpha = alpha)
            plt.plot(df['t'], df[col], s , label = fr"{name}", alpha = alpha, markersize = markersize, color = color)


open_file("wynik/euler_result_2y.txt", alpha=0.9, color = "red", name = "Metoda Eulera")
open_file("wynik/mid_point_result_2y.txt", alpha=0.9, color = "green", name = "Zmodyfikowana metoda Eulera", markersize=10)
open_file("wynik/heun_result_2y.txt", alpha=0.9, color = "blue", name = "Matoda Heuna")
open_file("wynik/rk4_result_2y.txt", alpha=0.5, color = "magenta", name = "RK4", s = "v", markersize=8)

# open_file("wynik/fehlberg_result_2y.txt", alpha=1, name = "RKF", color="green")
# open_file("wynik/cash_karp_result_2y.txt", alpha=1, name = "RKF-CK", color = (0/255, 87/255, 59/255, 1))

open_file("wynik/solution_2y.txt", s = "-", name = "Rozwiązanie analityczne")
open_sol("wynik/solution_sin.txt", s = "-", name = "Rozwiązanie analityczne")

# open_file("wynik/euler_result_sin.txt", alpha=0.8, color = "red", name = "Metoda Eulera")
# open_file("wynik/mid_point_result_sin.txt", alpha=0.8, color = "green", name = "Zmodyfikowana metoda Eulera")
# open_file("wynik/heun_result_sin.txt", alpha=0.8, name = "Metoda Heuna", color = "blue")
# open_file("wynik/rk4_result_sin.txt",  alpha=1, name = "RK4", color="magenta", s = "v", markersize = 8)

# open_file("wynik/cash_karp_result_sin.txt", alpha=0.5, name = "stare")
# open_file("wynik/fehlberg_result_y'''.txt", s = "s" ,alpha=1, name = "RKF")
# open_file("wynik/cash_karp_result_y'''.txt", s = "s" ,alpha=1, name = "RKF-CK", color = (0/255, 87/255, 59/255, 1))



# open_file("wynik/euler_result_bez.txt", alpha=0.8, color = "red", name = "Metoda Eulera")
# open_file("wynik/mid_point_result_bez.txt", alpha=0.8, color = "green", name = "Zmodyfikowana metoda Eulera")
# open_file("wynik/heun_result_bez.txt", alpha=0.8, name = "Metoda Heuna", color = "blue")
# open_file("wynik/rk4_result_bez.txt",  alpha=1, name = "RK4", color="magenta", s = "v", markersize = 8)

# open_file("wynik/fehlberg_result_bez.txt", s = "s" ,alpha=1, name = "RKF")
# open_file("wynik/cash_karp_result_bez.txt", s = "s" ,alpha=1, name = "RKF-CK", color = (0/255, 87/255, 59/255, 1))








# plt.xlim(9.5, 10.2)
plt.xlabel("t")
plt.ylabel("y")
plt.legend(loc='best', fontsize=20)
plt.grid(True)
# plt.tight_layout()
plt.show()
