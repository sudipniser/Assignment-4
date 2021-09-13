'''
#A function to read the data from the files to obtain the matrix
#The file must be created in the following format, 
    #only one row must be written in a line, with elements separated by commas,
    #no comma is allowed to be placed at the end of a line
    #[1 0 0] is entered as 1,0,0
     [0 1 0]               0,1,0
     [0 0 1]               0,0,1
'''
def rd_mtrx(file_name):
    M=[]
    with open(str(file_name),'r') as file:
        for line in file:
            M+=[[]]
            fields= line.split(',')
            for num in fields:
                M[-1]+=[float(num)]
    return(M)

'''
#Creating a matrix class that will facilitate the handling for matrices for any dimension
The matrix class will take input in list of lists format, e.g., the matrix 
[1 0 0] will be entered as [[1,0,0],[0,1,0],[0,0,1]]
[0 1 0]
[0 0 1]
'''
class Matrix:
    def __init__(self,list_of_lists):
        for i in range(len(list_of_lists)):
            if len(list_of_lists[i]) != len(list_of_lists[0]): #It must be a matrix afterall
                print("undefined matrix")
        self.lol=list_of_lists
        self.row_nos=len(list_of_lists)
        self.col_nos=len(list_of_lists[0])
    def __getitem__(self,indices):
        return(self.lol[indices[0]][indices[1]])
    def __setitem__(self,indices,newval):
        self.lol[indices[0]][indices[1]]=newval
    def subm(self,r1,r2,c1,c2):#method to find submatrix from row r1 to r2 anc column c1 to c2
        X=[]
        for i in range(r2-r1+1):
            X+=[self.lol[r1:r2+1][i][c1:c2+1]]
        return(Matrix(X))
    def show(self):
        for i in self.lol:
            print(i)
    def Transpose(self):#Method defining the Transpose of a matrix
        C=[]
        for _ in range(len(self.lol[0])):#Inverting the number of rows and columns
            C+=[[]]
        for i in range(len(self.lol)):
            for j in range(len(self.lol[i])):#appending the rows as columns
                C[j].append(self.lol[i][j])
        return(Matrix(C))


def Mat_product(X,Y,show=True):
    if isinstance(X,Matrix) == True and isinstance(Y,Matrix) ==True:
#Just a check to see if the matrices can be actually be multiplied
        if len(X.lol[0])!=len(Y.lol):
            print('Product not defined')
        else:
            m=len(X.lol)
            n=len(Y.lol[0])
            Z=[]
            for _ in range(m):#creating the new matrix
                Z+=[[]]
            for i in range(len(X.lol)):#populating the entries of the new matrix
                temp=0
#The next three lines multiply the necessary elements to find the i,l th element of the new matrix
                for l in range(len(Y.lol[0])):
                    for j in range(len(X.lol[i])):
                        temp+=X.lol[i][j]*Y.lol[j][l]
                    Z[i].append(temp)
                    temp=0
            if show==True:
                Matrix(Z).show()
            return(Matrix(Z))
    else:
        print("Not a matrix")

'''
A function to multiply a factor 'fact' to the Rth row(i.e., perform the operation: R-->R*fact)
'''
def multRow(mtrx,R,fact):
    for i in range(mtrx.col_nos):
        mtrx[R,i]=mtrx[R,i]*fact

'''
A function to multiply a factor 'fact' to each element of the R2th row and replace each element of R1th row 
with the sum of both.(i.e., perform the operation: R1-->R1+R2*fact)
'''

def rowAdd(mtrx,R1,R2,fact=1):
    for i in range(mtrx.col_nos):
        mtrx[R1,i]+=fact*mtrx[R2,i]

'''
This function carries out the partial pivot operation, it takes the list of lists as input and carries out 
partial pivot about the pivot number supplied.
'''

def partialPivot(mtrx,pivnum):
    piv=mtrx[pivnum,pivnum]
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    swp_nos=0
    if piv!=0:
        return(swp_nos)
    else:
        swp_nos+=1
        swap_row_ind=pivnum
        for i in range(pivnum,m):
            if mtrx[i,pivnum]>mtrx[swap_row_ind,pivnum]:
                swap_row_ind=i
        if mtrx[swap_row_ind,pivnum]==0:
            return("Unswappable zero pivot reached!")
        mtrx[pivnum,:],mtrx[swap_row_ind,:]=mtrx[swap_row_ind,:],mtrx[pivnum,:]
        return(swp_nos)

'''
Defining a function that will take an augmented matrix object, perform the Gauss jordan and return the 
Row reduced echelon matrix
'''

def GaussJordan(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    for i in range(m):
        pP=partialPivot(mtrx,i)
        if pP=="Unswappable zero pivot reached!":
            return("Unswappable zero pivot reached!")
        multRow(mtrx,i,1/mtrx[i,i])
        for j in range(m):
            if j!=i:
                rowAdd(mtrx,j,i,-mtrx[j,i])

'''
A function that takes in the Augmented matrix and returs the inverse.
Note: It does not return the augmented matrix but only the inverse
'''

def inverse(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    det_check=GaussJordan(mtrx)
    if det_check=="Unswappable zero pivot reached!":
        return('Determinant does not exist')
    else:
        C=[]
        for i in range(m):
            C+=[mtrx[i,m:]]
        return(Matrix(C))

'''
A function that calculates the determinant of a square matrix
'''

def Determinant(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    swp_nos=0
    for i in range(m):
        pP=partialPivot(mtrx,i)
        if pP=="Unswappable zero pivot reached!":
            return(0)
        swp_nos+=pP
        for j in range(i+1,m):
            rowAdd(mtrx,j,i,-mtrx[j,i]/mtrx[i,i])
    val=1
    for pivnum in range(m):
        val*=mtrx[pivnum,pivnum]
    return((-1)**(swp_nos)*val)

'''
A function that rounds off the elements of a matrix to two decimal places
'''

def round_Mat(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    for i in range(m):
        for j in range(n):
            mtrx[i,j]=float(format(mtrx[i,j],'.2f'))
    return(mtrx)

'''
#everything before this is mostly from previous assignment with some minor improvements
LU decompostion(Doolittle/l_{ii}=1)
#Takes a matrix input and decomposes the input into a lower triangular matrix and an upper triangluar matrix
#The original matrix is LOST
'''
def LUdecomp_Doolittle(mtrx):
    m=mtrx.row_nos         #number of rows
    n=mtrx.col_nos         #number of columns
    for k in range(m):
        #This section will check if the decompostion is possible
        #And if required will use partial pivoting to attempt factorisation
        if Determinant(mtrx.subm(0,k,0,k))==0:
            pP=partialPivot(mtrx,k) 
            if pP=="Unswappable zero pivot reached!":
                return("LU Decomposition not possible")   
    #This section caclulates the value of the elements in L and U       
    for i in range(1,m):                 #choosing row number
        for j in range(m):               #choosing column number
            if i<=j:                     #Calculating u_{ij}
                temp=0
                for k in range(i):
                    temp+=mtrx[i,k]*mtrx[k,j]
                mtrx[i,j]=mtrx[i,j]-temp #replacing the original matrix with the calculated value
            elif i>j:                    #Calculating l_{ij}
                temp=0
                for k in range(j):
                    temp+=mtrx[i,k]*mtrx[k,j]
                mtrx[i,j]=(mtrx[i,j]-temp)/mtrx[j,j]


'''
LU decompostion(Crout/u_{jj}=1)
#Takes a matrix input and decomposes the input into a lower triangular matrix and an upper triangluar matrix
#This section operates just like the previous one with a minor difference in the formula
#The original matrix is LOST
'''
def LUdecomp_Crout(mtrx):
    m=mtrx.row_nos  #number of rows
    n=mtrx.col_nos  #number of columns
    for k in range(m):
        #This section will check if the decompostion is possible
        #And if required, will use partial pivoting to attempt factorisation
        if Determinant(mtrx.subm(0,k,0,k))==0:
            pP=partialPivot(mtrx,k) 
            if pP=="Unswappable zero pivot reached!":
                return("LU Decomposition not possible")
    #This section caclulates the value of the elements in L and U
    for j in range(1,m):                           #choosing column number
        for i in range(m):                         #choosing row number
            if i>=j:                               #Calculating u_{ij}
                temp=0
                for k in range(j):
                    temp+=mtrx[i,k]*mtrx[k,j]
                mtrx[i,j]=mtrx[i,j]-temp            #Notice the difference in formula
            if i<j:                                 #calculating l_{ij}
                temp=0
                for k in range(i):
                    temp+=mtrx[i,k]*mtrx[k,j]
                mtrx[i,j]=(mtrx[i,j]-temp)/mtrx[i,i]#Notice the difference in formula

'''
Cholesky decompostion
#Only works for Hermitian and positive definite matrices
#To compute square root, 'math' library has been imported
#Takes a matrix input and decomposes the input into a lower triangular matrix and an upper triangluar matrix
#This section operates just like the previous ones with a minor difference in the formula
#The original matrix is LOST
'''

def CholDecomp(mtrx):
    import math
    m=mtrx.row_nos         #number of rows
    n=mtrx.col_nos         #number of columns
    for k in range(m):
        #This section will check if the decompostion is possible
        #And if required, will use partial pivoting to attempt factorisation
        if Determinant(mtrx.subm(0,k,0,k))==0:
            pP=partialPivot(mtrx.lol,k) 
            if pP=="Unswappable zero pivot reached!":
                return("Cholesky Decomposition not possible") 
    #This section caclulates the value of the elements in L and U
    #we dont need to calculate every combination of i,j since the matrix is symmetric
    for i in range(m):
        for j in range(i,m):
            if i==j:                   #diagonal elements
                temp=0
                for k in range(i):
                    temp+=mtrx[i,k]**2
                mtrx[i,i]=math.sqrt(mtrx[i,i]-temp)
            if i<j:                    #calculating off diagonal elements
                temp=0
                for k in range(i):
                    temp+=mtrx[i,k]*mtrx[k,j]
                mtrx[i,j]=(mtrx[i,j]-temp)/mtrx[i,i]
                mtrx[j,i]=mtrx[i,j]    #since the matrix is symmetric


'''
backward and forward substitution for Doolitlle's method
'''

def backsubs_Doolittle(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    X=[]
    for _ in range(m):
        X+=[[]]
    for i in range(m-1,-1,-1):#moving backwards
        temp=0
        for j in range(i+1,m):
            temp+=mtrx[i,j]*X[j][0]
        X[i]=[(mtrx[i,-1]-temp)/mtrx[i,i]]
    return(Matrix(X))
def forsubs_Doolittle(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    X=[]
    for _ in range(m):
        X+=[[]]
    for i in range(m):#moving forward
        temp=0
        for j in range(i):
            temp+=mtrx[i,j]*X[j][0]
        X[i]=[(mtrx[i,-1]-temp)]
    return(Matrix(X))


'''
backward and forward substitution for Crout's method
'''

def backsubs_Crout(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    X=[]
    for _ in range(m):
        X+=[[]]
    for i in range(m-1,-1,-1):#moving backwards
        temp=0
        for j in range(i+1,m):
            temp+=mtrx[i,j]*X[j][0]
        X[i]=[(mtrx[i,-1]-temp)]
    return(Matrix(X))
def forsubs_Crout(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    X=[]
    for _ in range(m):
        X+=[[]]
    for i in range(m):#moving forward
        temp=0
        for j in range(i):
            temp+=mtrx[i,j]*X[j][0]
        X[i]=[(mtrx[i,-1]-temp)/mtrx[i,i]]
    return(Matrix(X))


'''
backward and forward substitution for Cholesky decomposition
'''
def backsubs_chol(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    X=[[] for _ in range(m)]
    for i in range(m-1,-1,-1):#moving backwards
        temp=0
        for j in range(i+1,m):
            temp+=mtrx[i,j]*X[j][0]
        X[i]=[(mtrx[i,-1]-temp)/mtrx[i,i]]
    return(Matrix(X))
def forsubs_chol(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    X=[[] for _ in range(m)]
    for i in range(m):
        temp=0
        for j in range(i):#moving forward
            temp+=mtrx[i,j]*X[j][0]
        X[i]=[(mtrx[i,-1]-temp)/mtrx[i,i]]
    return(Matrix(X))

'''
Solve a system of equation using LU decomposition(Doolittle method)
AX=b
#takes A and b(separately) as input and returns X
#Original matrices are preserved
'''

def solve_system_Doolittle(A,b):
    m=A.row_nos   #number of rows
    n=A.col_nos   #number of columns
    Aug=[]
    for i in range(m):#making the augmented matrix
        Aug.append(A[i,:]+b[i,:])
    LUb=Matrix(Aug)
    check=LUdecomp_Doolittle(LUb)#LU decompostion
    if check=="LU Decomposition not possible":
        return("No unique solution exists")
    y=forsubs_Doolittle(LUb)#forward substitution to solve Ly=b
    Uy=[]
    for i in range(m):
        Uy.append(LUb[i,:-1]+y[i,:])
    x=backsubs_Doolittle(Matrix(Uy))#backward substitution to solve UX=y
    return(x)
'''
Solve a system of equation using LU decomposition(Crout's method)
AX=b
#takes A and b(separately) as input and returns X
#Original matrices are preserved
'''

def solve_system_Crout(A,b):
    m=A.row_nos   #number of rows
    n=A.col_nos   #number of columns
    Aug=[]
    for i in range(m):
        Aug.append(A[i,:]+b[i,:])
    LUb=Matrix(Aug)#making the augmented matrix
    check=LUdecomp_Crout(LUb)
    if check=="LU Decomposition not possible":
        return("No unique solution exists")
    y=forsubs_Crout(LUb)#forward substitution to solve Ly=b
    Uy=[]
    for i in range(m):
        Uy.append(LUb[i,:]+y[i,:])
    x=backsubs_Crout(Matrix(Uy))#backward substitution to solve UX=y
    return(x)

'''
Solve a system of equation using Cholesky decomposition
AX=b
#A is hermitian and positive definite
#takes A and b(separately) as input and returns X
#Original matrices are preserved
'''

def solve_system_chol(A,b):
    m=A.row_nos   #number of rows
    n=A.col_nos   #number of columns
    Aug=[]
    for i in range(m):
        Aug.append(A[i,:]+b[i,:])
    LUb=Matrix(Aug)#making the augmented matrix
    check=CholDecomp(LUb)
    if check=="Cholesky Decomposition not possible":
        return("No unique solution exists")
    y=forsubs_chol(LUb)#forward substitution to solve Ly=b
    Uy=[]
    for i in range(m):
        Uy.append(LUb[i,:]+y[i,:])
    x=backsubs_chol(Matrix(Uy))#backward substitution to solve UX=y
    return(x)

'''
Inverse using LU decomposition
#Calculates inverse by solving AX=I for each row of I
#Takes only the matrix (not augmented) as input and returns the inverse
#original matrix may be LOST
'''

def LUinverse(mtrx):
    m=mtrx.row_nos   #number of rows
    n=mtrx.col_nos   #number of columns
    I=[]             #Indentity matrix
    inv=[]
    for i in range(m):
        temp=[float(j==i) for j in range(m) ]
        I+=[temp]
    for i in range(m):
        temp=solve_system_Doolittle(mtrx,Matrix([[ele[i]] for ele in I]))
        if temp=="No unique solution exists":
            return("inverse does not exist")
        inv+=temp.Transpose().lol
    return(Matrix(inv).Transpose())