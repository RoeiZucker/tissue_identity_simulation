from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(3)
values = {'Skin': (0.7469556156341755, 0.08832391849287316), 'Sperm': (0.8768583221892472, 0.06728471180982226), 'Urine': (0.9159405960251148, 0.05661383700013867)}
results = [0.746 * 100, 0.876 * 100, 0.913 * 100]
stderr = [0.088,0.877,0.057]



# {'Skin': 0.7460792517158649, 'Sperm': 0.8755033021820008, 'Urine': 0.9130585836618551}

def millions(x, pos):
    'The two args are the value and tick position'
    return '$%1.1fM' % (x * 1e-6)

formatter = FuncFormatter(millions)

fig, ax = plt.subplots()
# ax.yaxis.set_major_formatter(formatter)
ax.set_ylim([0, 100])
plt.title("Correct allocation rate per tissue")
ax.set_ylabel("Correct allocation percentage")

plt.xticks(x, ('Person A - Skin', 'Person B - Sperm', 'Person C - Urine'))
# ind = np.arange(len(results))  # the x locations for the groups
ind = np.arange(3)
rects = ax.bar(ind, results)
counter = 0
for rect in rects:
    height = rect.get_height()
    ax.annotate('{}'.format(height),
                xy=(rect.get_x() + rect.get_width() / 2 - 0.075, height + 2.5))
    plt.text(rect.get_x() + rect.get_width() / 2.0 + 0.24, height + 1.8, "±" + str(stderr[counter]), ha='center', va='bottom')
    counter+=1

plt.bar(x, results, color=['C2', 'C1', 'C0'])

plt.show()
