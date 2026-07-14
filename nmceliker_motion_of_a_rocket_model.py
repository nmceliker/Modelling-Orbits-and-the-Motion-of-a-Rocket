import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

# constants
G = 6.67430e-11  # gravitational constant in m^3 kg^-1 s^-2
M_E = 5.972e24  # mass of the earth in kg
M_M = 7.342e22  # mass of the moon in kg
r_E = 6371e3  # radius of the earth in m
r_M = 1737.4e3  # radius of the moon in m
distance_M = 384.4e6  # distance between the moon and the earth
m_rocket = 2420  # mass of the explorer in kg

# variables
vi = 10580  # initial velocity of the rocket in m/s   --> (11200m/s is the escape velocity)
xi = -7000e3  # initial x position in m
yi = 0  # initial y position in m
theta = np.pi / 2  # angle of launch in rad

# creating the time array
ti = 0  # initial time in s
tf = 2e6  # final time in s
dt = 10  # time step


# trajectory with the earth and the moon
def trajectory(t, y):
    x, y, vx, vy = y
    R_E = np.sqrt(x**2 + y**2)
    ax_E = -(G * M_E * x) / (R_E ** (3))
    ay_E = -(G * M_E * y) / (R_E ** (3))
    ax, ay = ax_E, ay_E
    if orbit == False:
        R_M = np.sqrt((x - distance_M) ** 2 + y**2)
        ax_M = -(G * M_M * (x - distance_M)) / (np.abs(R_M) ** (3))
        ay_M = -(G * M_M * y) / (np.abs(R_M) ** (3))
        ax += ax_M
        ay += ay_M
    return [vx, vy, ax, ay, ax_M, ay_M]


def crash_E(t, y):
    x, y, _, _ = y
    return x**2 + y**2 - r_E**2


def crash_M(t, y):
    x, y, _, _ = y
    return (x - distance_M) ** 2 + y**2 - r_M**2


def check_coor(x, y):
    valid_E = np.sqrt(x**2 + y**2)
    valid_M = np.sqrt((x - distance_M) ** 2 + y**2)
    if valid_E <= r_E:
        return False, True
    elif valid_M <= r_M:
        return True, False
    else:
        return True, True


def t_1st_orbit(solution):
    pos_y_index = np.where(solution.y[1] > 0)[0]
    if len(pos_y_index) > 0:
        close_x0 = pos_y_index[np.argmin(np.abs(solution.y[0][pos_y_index]))]

        x_index_sorted = np.argsort(np.abs(solution.y[0][pos_y_index]))
        close_x0_2 = pos_y_index[x_index_sorted[1]]

        t1_close = solution.t[close_x0]
        t2_close = solution.t[close_x0_2]

        t_orbit = t2_close - t1_close
        return t_orbit
    else:
        return None


MyInput = "0"
while MyInput != "q":
    MyInput = input('Enter a choice, "a", "b" or "q" to quit: ')
    print("You entered the choice: ", MyInput)
    if MyInput == "a":
        print("You have chosen part (a): simulation of an orbit")
        orbit = True
        coor_yn = input("Would you like to specify x and y co-ordinates? (y/n): ")
        if coor_yn == "y":
            Input_coor = input(
                "Enter x, y values (floating point position in meters, e.g. -6500000, 0): "
            )
            Input_coor = Input_coor.replace(" ", "").split(",")
            try:
                xi, yi = float(Input_coor[0]), float(Input_coor[1])
            except:
                print("Invalid entry, missing co-ordinates or wrong format.")
                continue
            xi, yi = float(Input_coor[0]), float(Input_coor[1])
            valid_E, valid_M = check_coor(xi, yi)
            if not valid_E and valid_M:
                print("Invalid input, coordinates lie within the Earth.")
                continue
        elif coor_yn == "n":
            print(f"Preset values of {xi}, {yi} are going to be used.")
        else:
            print("Invalid choice.")
        vi_yn = input("Would you like to set an initial velocity in m/s (y/n): ")
        if vi_yn == "y":
            Input_vi = input("Enter a value for initial velocity, vi in m/s: ")
            vi = float(Input_vi)
        elif vi_yn == "n":
            print(
                f"You have chosen not to specify an initial velocity. Previously set vi of {vi} m/s is going to be used."
            )
        else:
            print("Invalid choice.")
        theta_yn = input("Would you like to specify a launch angle (y/n): ")
        if theta_yn == "y":
            theta_input = input("Enter a value for the launch angle in degrees: ")
            theta = float(float(theta_input) * (np.pi / 180))
        elif theta_yn == "n":
            print(
                "You have chosen not to specify a launch angle, previously set angle of 90 degrees is going to be used."
            )
        else:
            print("Invalid choice.")
        tf_yn = input("Would  you like to specify a simulation length (y/n): ")
        if tf_yn == "y":
            tf = float(input("Enter a simulation length in seconds: "))
        elif tf_yn == "n":
            print(
                f"You have chosen not to specify a simulation length, previously set value of {tf} s is going to be ued."
            )
        else:
            print("Invalid choice.")
        KE_yn = input(
            "Would you like to plot KE and PE against time plots for the explorer (y/n): "
        )
        if KE_yn == "y":
            KE_yn_TF = True
        elif KE_yn == "n":
            print("You have  chosen not to produce plots for KE and PE against time.")
        else:
            print("Invalid input.")

        # x and y components of the initial velocity in m/s
        vxi = vi * np.cos(theta)
        vyi = vi * np.sin(theta)
        initial_conditions = [xi, yi, vxi, vyi]

        t_total = (ti, tf)  # total time of the simulation seconds
        t_arr = np.arange(
            t_total[0], t_total[1] + dt, dt
        )  # an array for the tiem w/ dt

        # to solve the ODE using solve_ivp within the scipy library
        solution = scipy.integrate.solve_ivp(
            trajectory,
            t_total,
            initial_conditions,
            method="RK45",
            events=[crash_E, crash_M],
            rtol=1e-7,
            atol=1e-10,
            t_eval=t_arr,
        )

        possible_events = [set(zip(events, events)) for events in solution.t_events]

        if any(("crash_E", "crash_E") in events for events in possible_events):
            print("Rocket crashed into the Earth.")
        else:
            print("Rocket completed its trajectory without crashing.")

        # plots the trajectory
        plt.figure(figsize=(8, 6))
        plt.plot(solution.y[0], solution.y[1], label="Rocket")
        plt.plot(0, 0, "o", color="blue", markersize=15, label="Earth")
        plt.title("Trajectory of a Rocket Around the Earth")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.axis("equal")
        plt.grid(True)
        plt.legend()
        plt.show()

        t_orbit = t_1st_orbit(solution)
        t_orbit_day = round(t_orbit / 86400, 2)
        if t_orbit is not None:
            print(f"First orbit was completed in {t_orbit}s or {t_orbit_day} day(s).")
        else:
            print("No period for the first orbit was calculated.")

        if KE_yn_TF == True:
            KE_rocket = 0.5 * m_rocket * (solution.y[2] ** 2 + solution.y[3] ** 2)
            PE_rocket_E_total = m_rocket * (
                solution.y[2] * solution.y[0] + solution.y[3] * solution.y[1]
            )

            plt.figure(figsize=(10, 6))
            plt.plot(solution.t, KE_rocket, color="blue")
            plt.plot(solution.t, PE_rocket_E_total, color="red")
            plt.title("Kinetic Energy of the Rocket over Time")
            plt.xlabel("Time (seconds)")
            plt.ylabel("Kinetic Energy (Joules)")
            plt.grid(True)
            plt.show()

    elif MyInput == "b":
        print("You have chosen part (b): shooting the moon")
        orbit = False
        coor_yn = input("Would you like to specify x and y co-ordinates? (y/n): ")
        if coor_yn == "y":
            Input_coor = input(
                "Enter x, y values (floating point position in meters, e.g. -6500000, 0): "
            )
            Input_coor = Input_coor.replace(" ", "").split(",")
            try:
                xi, yi = float(Input_coor[0]), float(Input_coor[1])
            except:
                print("Invalid entry, missing co-ordinates or wrong format.")
                continue
            xi, yi = float(Input_coor[0]), float(Input_coor[1])
            valid_E, valid_M = check_coor(xi, yi)
            if not valid_E and valid_M:
                print("Invalid input, coordinates lie within the Earth.")
                continue
            elif valid_E and not valid_M:
                print("Invalid input, coordinates lie within the Moon.")
                continue
        elif coor_yn == "n":
            print(f"Preset values of {xi}, {yi} are going to be used.")
        else:
            print("Invalid choice.")
        vi_yn = input("Would you like to set an initial velocity in m/s (y/n): ")
        if vi_yn == "y":
            Input_vi = input("Enter a value for initial velocity, vi in m/s: ")
            vi = float(Input_vi)
        elif vi_yn == "n":
            print(
                f"You have chosen not to specify an initial velocity. Previously set vi of {vi} m/s is going to be used."
            )
        else:
            print("Invalid choice.")
        theta_yn = input("Would you like to specify a launch angle (y/n): ")
        if theta_yn == "y":
            theta_input = input("Enter a value for the launch angle in degrees: ")
            theta = float(float(theta_input) * (np.pi / 180))
        elif theta_yn == "n":
            print(
                "You have chosen not to specify a launch angle, previously set angle of 90 degrees is going to be used."
            )
        else:
            print("Invalid choice.")
        tf_yn = input("Would  you like to specify a simulation length (y/n): ")
        if tf_yn == "y":
            tf = float(input("Enter a simulation length in seconds: "))
        elif tf_yn == "n":
            print(
                f"You have chosen not to specify a simulation length, previously set value of {tf} s is going to be ued."
            )
        else:
            print("Invalid choice.")
        KE_yn = input(
            "Would you like to plot KE and PE against time plots for the explorer (y/n): "
        )
        if KE_yn == "n":
            print("You have  chosen not to produce plots for KE and PE against time.")
        # x and y components of the initial velocity in m/s
        vxi = vi * np.cos(theta)
        vyi = vi * np.sin(theta)
        initial_conditions = [xi, yi, vxi, vyi]

        t_total = (ti, tf)  # total time of the simulation seconds
        t_arr = np.arange(
            t_total[0], t_total[1] + dt, dt
        )  # an array for the tiem w/ dt

        # to solve the ODE using solve_ivp within the scipy library
        solution = scipy.integrate.solve_ivp(
            trajectory,
            t_total,
            initial_conditions,
            method="RK45",
            events=[crash_E, crash_M],
            rtol=1e-7,
            atol=1e-10,
            t_eval=t_arr,
        )

        possible_events = [set(zip(events, events)) for events in solution.t_events]

        if any(("crash_E", "crash_E") in events for events in possible_events):
            print("Rocket crashed into the Earth.")
        elif any(("crash_M", "crash_M") in events for events in possible_events):
            print("Rocket crashed into the Moon.")
        else:
            print("Rocket completed its trajectory without crashing.")

        # plots the trajectory
        plt.figure(figsize=(8, 6))
        plt.plot(solution.y[0], solution.y[1], label="Rocket")
        plt.plot(0, 0, "o", color="blue", markersize=10, label="Earth")
        plt.plot(distance_M, 0, "o", color="gray", markersize=5, label="Moon")
        plt.title("Trajectory of a Rocket Launched to the Moon")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.axis("equal")
        plt.grid(True)
        plt.legend()
        plt.show()

        d_to_M = np.sqrt((solution.y[0] - distance_M) ** 2 + solution.y[1] ** 2)
        if d_to_M.size > 0:
            d_to_M_index = np.argmin(d_to_M)
            closest_to_M = d_to_M[d_to_M_index]
            print(f"Closest the rocket got to the moon was {closest_to_M}.")
        else:
            print("Distance to the Moon was not calculated.")

        t_orbit = t_1st_orbit(solution)
        t_orbit_day = round(t_orbit / 86400, 2)
        if t_orbit is not None:
            print(f"First orbit was completed in {t_orbit}s or {t_orbit_day} day(s).")
        else:
            print("No period for the first orbit was calculated.")

        if KE_yn_TF == True:
            KE_rocket = 0.5 * m_rocket * (solution.y[2] ** 2 + solution.y[3] ** 2)
            PE_rocket_E_E = m_rocket * (
                solution.y[2] * solution.y[0] + solution.y[3] * solution.y[1]
            )
            PE_rocket_E_M = m_rocket * (
                solution.y[4] * (solution.y[0] - distance_M)
                + solution.y[5] * solution.y[1]
            )
            PE_rocket_tot = PE_rocket_E_E + PE_rocket_E_M

            plt.figure(figsize=(10, 6))
            plt.plot(solution.t, KE_rocket, color="blue")
            plt.plot(solution.t, PE_rocket_tot, color="red")
            plt.title("Kinetic Energy of the Rocket over Time")
            plt.xlabel("Time (seconds)")
            plt.ylabel("Kinetic Energy (Joules)")
            plt.grid(True)
            plt.show()

    elif MyInput != "q":
        print("This is not a valid choice")
print("You have chosen to finish - goodbye.")
