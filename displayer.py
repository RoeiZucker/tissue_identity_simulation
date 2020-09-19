import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

corr = {'Urine': 0, 'Sperm': 0, 'Skin': 15.27, 'Saliva_edited': 14.23, 'Whole_blood': 12.42}
incorr = {'Urine': 2.76, 'Sperm': 3.65, 'Skin': 6.32, 'Saliva_edited': 5.06, 'Whole_blood': 1.67}
labels = ['Urine', 'Sperm', 'Skin', 'Saliva', 'Whole Blood']
correct_means = [corr["Urine"], corr["Sperm"], corr["Skin"], corr["Saliva_edited"], corr["Whole_blood"]]
incorrect_means = [incorr["Urine"], incorr["Sperm"], incorr["Skin"], incorr["Saliva_edited"], incorr["Whole_blood"]]

x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width / 2, correct_means, width, label='Correctly guessed')
rects2 = ax.bar(x + width / 2, incorrect_means, width, label='Incorrectly guessed')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('mean number of discoveries')
ax.set_title('Assigned tissues for 100 runs for 1/3 Skin, Saliva and Whole Blood')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)

fig.tight_layout()

plt.show()
