#!!!!!!!!!!!!!!#
#已经成功运行到Sn=7！！！！！
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!注意每个（in）矩阵除了（n-1，n）外都未验证


'''
由已知的杨盘，杨图，分支律产生Sn的交换矩阵（ij）（in）的程序
'''


#输入Sn in [3,10]，外部输入到Sn=2
Sn=7

#生成其他顶层变量
print('计算Sn='+str(Sn)+'的（ij）（in）交换矩阵程序现在开始！')
St=Sn-1

#import索引
#os；pickle；math；numpy
import pickle,os,math,copy
import numpy as np
#import scipy as sp

#生成Sn对应的(ij),(in)列表，目录，文件夹
Allij=[[Sn-i-1,Sn-i]for i in range(Sn-1)]
IJpath=['C:\\CG\\YamanouchiMatrix(ij)\\S'+str(Sn)+'\\ij'+str(i)+'\\' for i in Allij]
Allin=[[Sn-i-1,Sn]for i in range(Sn-1)]
INpath=['C:\\CG\\YamanouchiMatrix(in)\\S'+str(Sn)+'\\in'+str(i)+'\\' for i in Allin]
print(IJpath,INpath)
for i in IJpath:
    os.makedirs(i)
for i in INpath:
    os.makedirs(i)
print('Sn='+str(Sn)+'的（ij）（in）目录生成结束')

#读取存在C:\CG\r\Sn.pkl中的杨盘[ri]   #这里是Sn下所有的【ri】,所以[ri]是双列表哦！
Myfile=open('C:\\CG\\r\\S'+str(Sn)+'.pkl','rb')
r=pickle.load(Myfile)
Myfile.close()
print('读取[ri]成功，让我们欢快地开始造Sn=S'+str(Sn)+'的Yamanouchi矩阵吧！')

#定义函数
#寻找矩阵c中元素d坐标的函数
def FindSeat(Matrix,value):
    '''
    寻找矩阵c中元素d坐标的函数
    :param c:Mateix
    :param d: value
    :return: Seat of value
    '''
    x=0
    y=0
    for i in Matrix:
        try:
            y=i.index(value)
        except:
            x+=1
        else:
            break
    Seat=[x,y]
    return Seat
    #end
#Function1算当前元素自己的交换矩阵元
def Function1(YangPic):
    SeatSt=FindSeat(YangPic,St)
    SeatSn=FindSeat(YangPic,Sn)
    if SeatSt[0]==SeatSn[0]:
        return 1
    elif SeatSt[1]==SeatSn[1]:
        return -1
    else:
        Xigema=SeatSt[0]-SeatSt[1]+SeatSn[1]-SeatSn[0]
        Yuansu=1/Xigema
        return Yuansu
    #end
#Function2算YangPic经过（n-1，n）变换后的对应杨图的相对位置
def Function2(YangPic,r):
    Table=copy.deepcopy(YangPic)
    SeatSt=FindSeat(YangPic, St)
    SeatSn=FindSeat(YangPic, Sn)
    print('Table=',Table,'SeatSt=',SeatSt,'SeatSn=',SeatSn)
    Table[SeatSn[0]][SeatSn[1]]=St
    Table[SeatSt[0]][SeatSt[1]]=Sn
    print('Table后=', Table, 'SeatSt=', SeatSt, 'SeatSn=', SeatSn)
    number=0
    dissR=copy.deepcopy(r)
    for sn in range(Sn,0,-1):
        if sn==1:
            number+=1
        else:
            Seatsn=FindSeat(Table,sn)
            Myfile=open('C:\\CG\\BranchLaw\\S'+str(sn)+'\\BeforeS\\'+str(dissR)+'.pkl','rb')
            BeforeR=pickle.load(Myfile)
            Myfile.close()
            print('Seatsn=',Seatsn)
            if dissR[Seatsn[0]]==1:
                dissR=dissR[:-1]
            else:
                dissR[Seatsn[0]]+=-1#一定要在loadBeforeR之后再改为删去sn的杨盘
            k=BeforeR.index(dissR)
            if k>0:
                for k in range(k):
                    Myfile=open('C:\\CG\\YangPic\\S'+str(sn-1)+'\\'+str(BeforeR[k])+'\\Mm.txt','r')
                    M=Myfile.readline()
                    Myfile.close()
                    M=int(M)
                    number+=M
    return number
    #end


#正式开始遍历[ri]造(ij)矩阵和(in)矩阵
print('请等待几分钟...')
for r in r:#下面的r就是列表了
    #读k和Before
    BLpath='C:\\CG\\BranchLaw\\S'
    Myfile=open(BLpath+str(Sn)+'\\'+'k\\'+str(r)+'.txt','r')
    k=Myfile.readline()
    Myfile.close()
    k=int(k)
    Myfile=open(BLpath+str(Sn)+'\\'+'BeforeS\\'+str(r)+'.pkl','rb')
    BeforeS=pickle.load(Myfile)
    Myfile.close()
    #读建立矩阵的维数，用以校验
    Yangpath='C:\\CG\\YangPic\\S'
    Myfile=open(Yangpath+str(Sn)+'\\'+str(r)+'\\Mm.txt','r')
    DR=Myfile.readline()
    Myfile.close()
    DR=int(DR)
    #建立DR维标准全零矩阵
    MatrixC=[[0 for i in range(DR)]for j in range(DR)]

    #开始做(n-1,n)
    IN=[Sn-1,Sn]
    MatrixIN=copy.deepcopy(MatrixC)#他存活到被下一个in（n-2，n）刷新，过程中要与ij（n-2,n-1）做一次乘法
    m=0   #过程中表示每个分快开始位置
    for i in range(k):
        MyfileMm=open('C:\\CG\\YangPic\\S'+str(St)+'\\'+str(BeforeS[i])+'\\Mm.txt','r')
        M=MyfileMm.readline()
        MyfileMm.close()
        M=int(M)
        print('M=',M,'M路径','C:\\CG\\YangPic\\S'+str(St)+'\\'+str(BeforeS[i])+'\\Mm')
        for j in range(M):
            Myfile=open(Yangpath+str(Sn)+'\\'+str(r)+'\\'+str(j+m+1)+'.pkl','rb')
            YangPic=pickle.load(Myfile)
            Myfile.close()
            print('YangPic=',YangPic,'j',j,'YangPic路径',Yangpath+str(Sn)+'\\'+str(r)+'\\'+str(j+m+1)+'.pkl')
            MatrixIN[j+m][j+m]=Function1(YangPic)#算YangPic自己对应的矩阵元
            print('YangPic=',YangPic,'这个应该和上面一样才对')
            print('MatrixIN=',MatrixIN)
            if MatrixIN[j+m][j+m]!=-1 and MatrixIN[j+m][j+m]!=1:
                print('MatrixIN[j+m][j+m]='+str(MatrixIN[j+m][j+m])+'通过非1测试')
                SeatY=Function2(YangPic,r)#算YangPic经过（n-1，n）变换后的对应杨图的分块外位置【数学列坐标】
                print('YangPic=', YangPic, '这个应该和上面一样才对')
                SeatY+=-1#SeatY是python坐标【python列坐标】
                print('SeatY=',SeatY)

                MatrixIN[j+m][SeatY]=math.sqrt(1-MatrixIN[j+m][j+m]**2)#sqrt(x**2-1)/|x|=sqrt(1-1/x**2)
                print('MatrixIN=', MatrixIN)
        m+=M#更新分块
    if m==DR:  # 这个判断是一个分块和等于总数的安全性检查  ## 注意：前面仅仅计算了(n-1,n)这一个对换矩阵，它既是ij又是in
        #保存矩阵到（ij）
        Myfile=open(IJpath[0]+str(r)+'.pkl','wb')#pkl格式
        pickle.dump(MatrixIN,Myfile)
        Myfile.close()
        if Sn<6:#Sn小等于5的保存成txt副本便于查阅
            Myfile=open(IJpath[0]+str(r)+'.txt','w')
            Myfile.write(str(MatrixIN))
            Myfile.close()
        # 保存矩阵到（in）
        Myfile=open(INpath[0]+str(r)+'.pkl', 'wb')
        pickle.dump(MatrixIN,Myfile)
        Myfile.close()
        if Sn<6:#Sn小等于5的保存成txt副本便于查阅
            Myfile=open(INpath[0]+str(r)+'.txt', 'w')
            Myfile.write(str(MatrixIN))
            Myfile.close()
        #转换为mat格式保存和方便下一步计算
        MatrixIN=np.mat(MatrixIN)  # npy格式
        np.save(IJpath[0] + str(r) + '.npy', MatrixIN)
        np.save(INpath[0] + str(r) + '.npy', MatrixIN)

    #开始做剩下的(ij),(in)
    L=len(Allij)
    if L>1:
        for ij in Allij[1:]:  # 排除掉(n-1,n)外的剩余(ij)和(in)，Allij本省就是倒序哦！
            MatrixIJ=copy.deepcopy(MatrixC)
            m=0
            for i in range(k):
                Myfile=open('C:\\CG\\YamanouchiMatrix(ij)\\S'+str(St)+'\\ij'+str(ij)+'\\'+str(BeforeS[i])+'.pkl','rb')
                MatrixPart=pickle.load(Myfile)
                Myfile.close()
                MyfileMm=open('C:\\CG\\YangPic\\S'+str(St)+'\\'+str(BeforeS[i])+'\\Mm.txt','r')
                M=MyfileMm.readline()
                MyfileMm.close()
                M=int(M)
                for aa in range(M):  # 这个地方，就是分支律的再次应用，整块儿移动，也许可以用pandas再优化
                    for bb in range(M):
                        MatrixIJ[m+aa][m+bb]=MatrixPart[aa][bb]
                m+=M

            #保存（ij）到pkl和txt
            if m==DR:  # 这里，（ij）就是纯粹的使用Sn-1，by分支律，按块儿原封不动拼图的，一定要再确认一下论文，免得搞错
                #保存矩阵到（ij）
                Myfile=open('C:\\CG\\YamanouchiMatrix(ij)\\S'+str(Sn)+'\\ij'+str(ij)+'\\'+str(r)+'.pkl','wb')  # pkl格式
                pickle.dump(MatrixIJ, Myfile)
                Myfile.close()
                if Sn<6:#Sn小等于5的保存成txt副本便于查阅
                    Myfile=open('C:\\CG\\YamanouchiMatrix(ij)\\S'+str(Sn)+'\\ij'+str(ij)+'\\'+str(r)+'.txt','w')
                    Myfile.write(str(MatrixIJ))
                    Myfile.close()
                #转化（ij）
                MatrixIJ=np.mat(MatrixIJ)#npy格式
                np.save('C:\\CG\\YamanouchiMatrix(ij)\\S'+str(Sn)+'\\ij'+str(ij)+'\\'+str(r)+'.npy',MatrixIJ)

                #计算(in),保存为npy
                a=np.dot(MatrixIN,MatrixIJ)
                b=np.dot(a,MatrixIN)
                MatrixIN=b.copy()#np下的矩阵不支持[:]分片，要用copy（）函数
                IN=[ij[0],Sn]
                np.save('C:\\CG\\YamanouchiMatrix(in)\\S'+str(Sn)+'\\in'+str(IN)+'\\'+str(r)+'.npy',MatrixIN)

                #将b转化为list保存pkl和txt
                b=np.matrix.tolist(b)
                Myfile=open('C:\\CG\\YamanouchiMatrix(in)\\S'+str(Sn)+'\\in'+str(IN)+'\\'+str(r)+'.pkl','wb')#pkl格式
                pickle.dump(b, Myfile)
                Myfile.close()
                #
                if Sn<5:#Sn小等于4的保存成txt副本便于查阅
                    Myfile=open('C:\\CG\\YamanouchiMatrix(in)\\S'+str(Sn)+'\\in'+str(IN)+'\\'+str(r)+'.txt','w')
                    Myfile.write(str(b))
                    Myfile.close()

#### 整体评价，对于Sn，（Sn-1，Sn）是一个矩阵元一个矩阵元赋值的，其他（ij）是用分支律直接拼图的，其他（in）是矩阵乘法得来的

#打完收工！
print('打完收工，目录在readme')