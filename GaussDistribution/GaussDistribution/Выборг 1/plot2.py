import matplotlib.pyplot as plt
import csv
import seaborn as sns

with open('data_pure.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    data_pure = [float(row[0]) for row in spamreader]

with open('data_sym_dirt.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    data_asym_dirt = [float(row[0]) for row in spamreader]

with open('data_sym_clog.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    data_asym_clog = [float(row[0]) for row in spamreader]


sns.displot(data={'Чистая': data_pure, 'Засоренная': data_asym_dirt, 'Засоряющая': data_asym_clog}, kind = 'kde')
plt.show()