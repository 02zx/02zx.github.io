import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from torchvision import datasets
from torchvision.transforms import ToTensor
import numpy as np
from torch import optim
from sklearn.preprocessing import MinMaxScaler
torch.set_default_tensor_type(torch.DoubleTensor)
import matplotlib.pyplot as plt
import matplotlib as mpl


dire='/home/XZ/lammps/ML/ice/PH/ice7/100/4.8/425.6/video' #'./train_7comp' #10.5
name='bop' #'train'
data=np.loadtxt('{}/{}.dat'.format(dire,name))[:12288]
nfeatures=np.shape(np.array(data))[1]
batch_size = 32
validation_split = .3
shuffle_dataset = True
random_seed= 42


dataset_size = len(data)
indices = list(range(dataset_size))
split = int(np.floor(validation_split * dataset_size))
if shuffle_dataset :
    np.random.seed(random_seed)
    np.random.shuffle(indices)
train_indices, val_indices = indices[split:], indices[:split]

train_sampler = torch.utils.data.SubsetRandomSampler(train_indices)
valid_sampler = torch.utils.data.SubsetRandomSampler(val_indices)

train_loader = torch.utils.data.DataLoader(data, batch_size=batch_size, sampler=train_sampler)
validation_loader = torch.utils.data.DataLoader(data, batch_size=batch_size,sampler=valid_sampler)

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device} device")

class AE(nn.Module):

    def __init__(self, c):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(nfeatures, 10*nfeatures),
            nn.Tanh()
        )

        self.bottleneck = nn.Sequential(
            nn.Linear(10*nfeatures, c)
        )

        self.decoder = nn.Sequential(
            nn.Linear(c, 10*nfeatures),
            nn.Tanh(),
            nn.Linear(10*nfeatures, nfeatures)
        )

    def forward(self, x):
        x1 = self.encoder(x)
        x2 = self.bottleneck(x1)
        x3 = self.decoder(x2)
        return x3



all_loss = []
c_loss = []
count = 0
minerr = 10000


for c in range(1, 8):
    model = AE(c)
    lr = 0.001
    weight_decay = 1e-5
    loss = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    epochs = 80
    for i in range(epochs):
        for idx, batch in enumerate(train_loader):
            output = model(batch[0].view(-1, nfeatures))
            loss_value = loss(output, batch[0].view(-1, nfeatures))

            model.zero_grad()
            loss_value.backward()

            all_loss.append(loss_value)

            optimizer.step()

        # print(i)
        # print(f"All loss  {all_loss[-1]}")
        
    error = (sum(all_loss[-nfeatures+1:]) / (nfeatures-1)).detach().numpy()
    if error < minerr:
        minerr = error
        torch.save(model, 'm_{}'.format(count))
        torch.save(model.state_dict(), 'p_{}'.format(count))
        count = count + 1
        
    print(all_loss[1-nfeatures:])
    c_loss.append((sum(all_loss[1-nfeatures:]) / (nfeatures-1)).detach().numpy())



plt.figure(figsize=(4,3),dpi=300)
plt.plot(range(1, 8), c_loss,'-ob')
plt.xlabel('The number of bottleneck nodes', fontdict={'family':'Arial', 'weight': 'normal', 'size': 15})
plt.ylabel('Mean squared error', fontdict={'family':'Arial', 'weight': 'normal', 'size': 15})


best_c = 3
print(best_c)
model = torch.load('m{}'.format(best_c-1))
traindata = torch.tensor(data[train_indices])
valdata = torch.tensor(data[val_indices])
outputs = model(valdata)
mse = nn.MSELoss()
val_error = mse(valdata,outputs)
print(val_error)


def gaussian_noise(x,per):
    mu=0.0
    std = per * np.std(x) # for %10 Gaussian noise
    noise = np.random.normal(mu, std, size = x.shape)
    x_noisy = x + noise
    return x_noisy

train=data[train_indices]
train2=train.copy()
train2=torch.tensor(train2)
orig=model(train2)


error = []
for i in range(nfeatures):
  train = data[train_indices]
  train[:,i] = gaussian_noise(train[:,i],0.1)
  train = torch.tensor(train)
  outputs = model(train)
  mse = nn.MSELoss()
  error.append(mse(outputs,orig).detach().numpy())

#Perturbation
RI0 = []
for i in range(nfeatures):
  total = sum(error)
  RI0.append(error[i]/total)

error = []
for i in range(nfeatures):
  train = data[train_indices]
  train[:,i] = gaussian_noise(train[:,i],0.5)
  train = torch.tensor(train)
  outputs = model(train)
  mse = nn.MSELoss()
  error.append(mse(outputs,orig).detach().numpy())


RI1 = []
for i in range(nfeatures):
  total = sum(error)
  RI1.append(error[i]/total)

error = []
for i in range(nfeatures):
  train = data[train_indices]
  s = np.mean(train[:,i])
  train[:,i] = np.repeat(s,len(train[:,i]))
  train = torch.tensor(train)
  outputs = model(train)
  mse = nn.MSELoss()
  error.append(mse(outputs,orig).detach().numpy())

#improved stepwise
RI2 = []
for i in range(nfeatures):
  total = sum(error)
  RI2.append(error[i]/total)


BOPs = ["q2","q4","q6","q8","q10","q12"]
x_axis = np.arange(len(BOPs))
c1=(42/256,87/256,141/256,1/3)
c2=(42/256,87/256,141/256,2/3)
c3=(42/256,87/256,141/256,1)
plt.figure(figsize=(6,3),dpi=300)
ax = plt.subplot(111)
ax.bar(x_axis-0.2, RI0, width=0.2, label = 'Perturbation 10%', color='black', align='center')
ax.bar(x_axis, RI1, width=0.2, label = 'Perturbation 50%', color='grey', align='center')
ax.bar(x_axis+0.2, RI2, width=0.2, color='navy', align='center',label='Improved Stepwise')
plt.xticks(x_axis, BOPs,size=15)
plt.legend()
plt.ylabel('Relative Importance', fontdict={'family':'Arial', 'weight': 'normal', 'size': 15})


def get_features(name):
    def hook(model, input, output):
        features[name] = output.detach()
    return hook

model.bottleneck.register_forward_hook(get_features('bottleneck'))

data=np.loadtxt('{}/{}.dat'.format(dire,name))
batch_size = 32
print(data.shape)
train_loader = torch.utils.data.DataLoader(data, batch_size=batch_size)

PREDS = []
FEATS = []
# placeholder for batch features
features = {}

# loop through batches
for idx, inputs in enumerate(train_loader):
    # move to device
    inputs = inputs.to("cpu")

    # forward pass [with feature extraction]
    preds = model(inputs)

    # add feats and preds to lists
    PREDS.append(preds.detach().cpu().numpy())
    FEATS.append(features['bottleneck'].cpu().numpy())
PREDS = np.concatenate(PREDS)
FEATS = np.concatenate(FEATS)

print('- preds shape:', PREDS.shape)
print('- feats shape:', FEATS.shape)
nfeatures=FEATS.shape[1]


from sklearn.mixture import GaussianMixture
def GMM_c(model, X, n):
    prob = model.predict_proba(X)
    Sk = 0
    if n == 1:
        for k in range(len(prob)):
            if prob[k] != 0:
                Sk = Sk + prob[k] * np.log(prob[k])
            else:
                Sk = Sk

    else:
        for i in range(len(prob)):
            for j in range(n):
                if prob[i, j] != 0:
                    Sk = Sk + prob[i, j] * np.log(prob[i, j])
                else:
                    Sk = Sk

    # S = np.dot(prob,np.log(prob.T))
    # Sk = -sum(np.diagonal(S))
    return -Sk


n_components = np.arange(1, 10)
models = [GaussianMixture(n, covariance_type='full', random_state=0).fit(FEATS)
          for n in n_components]
BIC = [m.bic(FEATS) for m in models]
S = []
for n in np.arange(1,10):
    models = GaussianMixture(n_components=n, random_state=0,covariance_type='full').fit(FEATS)
    t = GMM_c(models,FEATS,n)
    S.append(t)
ICL = []
for i in np.arange(len(BIC)):
    ICL.append(BIC[i] + 2 * S[i])



fig,ax=plt.subplots()
fig.figsize=(8,6)
fig.dpi=300
ax.plot(n_components[:6],BIC[:6],'-ro',label='BIC')
ax.set_ylabel('BIC', fontdict={'family':'Arial', 'weight': 'normal', 'size': 15})
plt.legend(prop={'size': 15}, loc="upper left")

plt.figure(figsize=(4,3),dpi=300)
plt.plot(n_components[:6],S[:6],'--bo',label='S')
plt.ylabel('Entropy', fontdict={'family':'Arial', 'weight': 'normal', 'size': 15})
plt.xlabel('Number of components', fontdict={'family':'Arial', 'weight': 'normal', 'size': 15})

BIC_min=3
GMM = GaussianMixture(n_components=BIC_min, random_state=0,covariance_type='full').fit(FEATS)
labels = GMM.predict(FEATS)
