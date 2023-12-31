plt.figure(figsize=(5,4),dpi=300)
r=np.arange(0,3.5,0.1)
s=4*np.pi*r**2
v=-(4*np.pi*r**3)/3
plt.gca().spines['right'].set_color('none')  
plt.gca().spines['bottom'].set_color('none')  
plt.gca().spines['top'].set_color('none')  
plt.axhline(y=0, color='black',linewidth=1)
plt.plot(r,s,'--k')
plt.plot(r,v,'-.k')
plt.plot(r,s+v,'-k')
plt.ylim(-20,20)
plt.xlim(0,3.5)
plt.xticks([])
plt.yticks([])
plt.legend(frameon=False)
plt.savefig('nuclei.png',format='png',bbox_inches='tight', transparent=True)
