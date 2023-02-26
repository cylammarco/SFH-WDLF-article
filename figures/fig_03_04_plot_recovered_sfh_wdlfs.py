import os

from matplotlib import pyplot as plt
import numpy as np
from spectres import spectres


# Load the GCNS data
gcns_wdlf = np.load(
    "pubgcnswdlf-h366pc-dpdf-samples-hp5-maglim80-vgen-grp-rdc-srt.npz"
)["data"]
n_bin_02 = 84
n_bin_05 = 32

h_gen_02, b_02 = np.histogram(
    gcns_wdlf["Mbol"],
    bins=n_bin_02,
    range=(1.9, 18.6),
    weights=0.01 / gcns_wdlf["Vgen"],
)
h_gen_05, b_05 = np.histogram(
    gcns_wdlf["Mbol"],
    bins=n_bin_05,
    range=(2.25, 18.25),
    weights=0.01 / gcns_wdlf["Vgen"],
)
e_gen_02, _ = np.histogram(
    gcns_wdlf["Mbol"],
    bins=n_bin_02,
    range=(1.9, 18.6),
    weights=0.01 / gcns_wdlf["Vgen"] ** 2.0,
)
e_gen_05, _ = np.histogram(
    gcns_wdlf["Mbol"],
    bins=n_bin_05,
    range=(2.25, 18.25),
    weights=0.01 / gcns_wdlf["Vgen"] ** 2.0,
)

bin_size_02 = b_02[1] - b_02[0]
bin_size_05 = b_05[1] - b_05[0]

mag_obs_02 = np.around(b_02[1:] - 0.5 * bin_size_02, 2)
obs_wdlf_02 = h_gen_02 / bin_size_02
obs_wdlf_err_02 = e_gen_02**0.5 / bin_size_02

mag_obs_05 = np.around(b_05[1:] - 0.5 * bin_size_05, 2)
obs_wdlf_05 = h_gen_05 / bin_size_05
obs_wdlf_err_05 = e_gen_05**0.5 / bin_size_05


# Load the pwdlfs
data = []
age_list_1 = np.arange(0.049, 0.100, 0.001)
age_list_2 = np.arange(0.100, 0.350, 0.005)
age_list_3 = np.arange(0.35, 15.01, 0.01)
age_list_3dp = np.concatenate((age_list_1, age_list_2))
age_list_2dp = age_list_3

age = np.concatenate((age_list_3dp, age_list_2dp))

for i in age_list_3dp:
    print(i)
    data.append(
        np.loadtxt(
            os.path.join(
                "output",
                f"montreal_co_da_20_K01_PARSECz0014_C08_{i:.3f}_Mbol.csv",
            ),
            delimiter=",",
        )
    )

for i in age_list_2dp:
    print(i)
    data.append(
        np.loadtxt(
            os.path.join(
                "output",
                f"montreal_co_da_20_K01_PARSECz0014_C08_{i:.2f}_Mbol.csv",
            ),
            delimiter=",",
        )
    )


mag_pwdlf = data[0][:, 0]

partial_age_02, solution_02 = np.load(
    "gcns_sfh_optimal_resolution_bin_02.npy"
).T
mag_obs_02, obs_wdlf_02, obs_wdlf_err_02 = np.load(
    "gcns_reconstructed_wdlf_optimal_resolution_bin_02.npy"
).T
partial_age_05, solution_05 = np.load(
    "gcns_sfh_optimal_resolution_bin_05.npy"
).T
mag_obs_05, obs_wdlf_05, obs_wdlf_err_05 = np.load(
    "gcns_reconstructed_wdlf_optimal_resolution_bin_05.npy"
).T
# Load the mapped pwdlf age-mag resolution
pwdlf_mapping_bin_02 = np.insert(np.load("pwdlf_bin_02_mapping.npy"), 0, 0)
pwdlf_mapping_bin_05 = np.insert(np.load("pwdlf_bin_05_mapping.npy"), 0, 0)

# Stack up the pwdlfs to the desired resolution
partial_wdlf_02 = []
partial_age_02 = []
for idx in np.sort(list(set(pwdlf_mapping_bin_02))):
    pwdlf_temp = np.zeros_like(mag_obs_02)
    age_temp = 0.0
    age_count = 0
    for i in np.where(pwdlf_mapping_bin_02 == idx)[0]:
        pwdlf_temp += spectres(mag_obs_02, mag_pwdlf, data[i][:, 1], fill=0.0)
        age_temp += age[i]
        age_count += 1
    partial_wdlf_02.append(pwdlf_temp)
    partial_age_02.append(age_temp / age_count)


partial_wdlf_05 = []
partial_age_05 = []
for idx in np.sort(list(set(pwdlf_mapping_bin_05))):
    pwdlf_temp = np.zeros_like(mag_obs_05)
    age_temp = 0.0
    age_count = 0
    for i in np.where(pwdlf_mapping_bin_05 == idx)[0]:
        pwdlf_temp += spectres(mag_obs_05, mag_pwdlf, data[i][:, 1], fill=0.0)
        age_temp += age[i]
        age_count += 1
    partial_wdlf_05.append(pwdlf_temp)
    partial_age_05.append(age_temp / age_count)


recomputed_wdlf_02 = np.nansum(
    solution_02 * np.array(partial_wdlf_02).T, axis=1
)
recomputed_wdlf_05 = np.nansum(
    solution_05 * np.array(partial_wdlf_05).T, axis=1
)


fig1, (ax1, ax_dummy1, ax2, ax_dummy2, ax3) = plt.subplots(
    nrows=5, ncols=1, figsize=(8, 10), height_ratios=(8, 1, 5, 1, 5)
)

ax_dummy1.axis("off")
ax_dummy2.axis("off")

ax1.errorbar(
    mag_obs_05,
    obs_wdlf_05,
    yerr=[obs_wdlf_err_05, obs_wdlf_err_05],
    fmt="+",
    markersize=5,
    label="Input WDLF",
)
ax1.plot(
    mag_obs_05,
    recomputed_wdlf_05
    / np.nansum(recomputed_wdlf_05)
    * np.nansum(obs_wdlf_05),
    label="Reconstructed WDLF",
)
ax1.set_xlabel(r"M${_\mathrm{bol}}$ / mag")
ax1.set_ylabel("log(arbitrary number density)")
ax1.set_xlim(0.0, 20.0)
ax1.set_ylim(1e-7, 5e-3)
ax1.set_yscale("log")
ax1.legend()
ax1.grid()

ax2.step(partial_age_05, solution_05, where="post")
ax2.grid()
ax2.set_xticks(np.arange(0, 15, 2))
ax2.set_xlim(0, 15)
ax2.set_ylim(bottom=0)
ax2.set_xlabel("Lookback time / Gyr")
ax2.set_ylabel("Relative Star Formation Rate")

ax3.step(partial_age_05, solution_05, where="post")
ax3.grid()
ax3.set_ylim(bottom=0)
ax3.set_xlabel("Lookback time / Gyr")
ax3.set_ylabel("Relative Star Formation Rate")
ax3.set_xscale("log")

plt.subplots_adjust(
    top=0.975, bottom=0.05, left=0.09, right=0.975, hspace=0.075
)

fig1.savefig("fig_03_gcns_reconstructed_wdlf_optimal_resolution_bin_05.png")


fig2, (ax1, ax_dummy1, ax2, ax_dummy2, ax3) = plt.subplots(
    nrows=5, ncols=1, figsize=(8, 10), height_ratios=(8, 1, 5, 1, 5)
)

ax_dummy1.axis("off")
ax_dummy2.axis("off")

ax1.errorbar(
    mag_obs_02,
    obs_wdlf_02,
    yerr=[obs_wdlf_err_02, obs_wdlf_err_02],
    fmt="+",
    markersize=5,
    label="Input WDLF",
)
ax1.plot(
    mag_obs_02,
    recomputed_wdlf_02
    / np.nansum(recomputed_wdlf_02)
    * np.nansum(obs_wdlf_02),
    label="Reconstructed WDLF",
)
ax1.set_xlabel(r"M${_\mathrm{bol}}$ / mag")
ax1.set_ylabel("log(arbitrary number density)")
ax1.set_xlim(0.0, 20.0)
ax1.set_ylim(1e-7, 5e-3)
ax1.set_yscale("log")
ax1.legend()
ax1.grid()

ax2.step(partial_age_02, solution_02, where="post")
ax2.grid()
ax2.set_xticks(np.arange(0, 15, 2))
ax2.set_xlim(0, 15)
ax2.set_ylim(bottom=0)
ax2.set_xlabel("Lookback time / Gyr")
ax2.set_ylabel("Relative Star Formation Rate")

ax3.step(partial_age_02, solution_02, where="post")
ax3.grid()
ax3.set_ylim(bottom=0)
ax3.set_xlabel("Lookback time / Gyr")
ax3.set_ylabel("Relative Star Formation Rate")
ax3.set_xscale("log")

plt.subplots_adjust(
    top=0.975, bottom=0.05, left=0.085, right=0.975, hspace=0.075
)

plt.savefig("fig_04_gcns_reconstructed_wdlf_optimal_resolution_bin_02.png")
