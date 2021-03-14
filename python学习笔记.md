1.  python中字符串和元组属于不可变序列（不能修改中间的指）

2.  python中允许多个变量指向同一个值，当其中一个变量修改值以后，其内存地址发生变化但不影响其他的变量

3.  python变量名区分大小写

4.  python里十六进制整数以0x头，八进制整数以0o开头，二进制整数以0b开头

5.  python的字符串格式化：

    ```python
    >>>a=3.6674
    >>>print('%7.3f'%a)
    '  3.667'
    >>>print('%06.2f'%a)
    '003.67' #小数点后四舍五入
    >>>print('%d:%c'%(65,65))
    '65:A'
    >>>print('my name is %s,and my age is %d'%('free'，19))
    'my name is free,and my age is 19'
    ```

6.  python要注意的运算符：

    | 运算符 |   功能   |
    | :----: | :------: |
    |   /    |  真除法  |
    |   //   |   整除   |
    |   *    | 乘，重复 |

7.  python内置函数：

    |         函数          |                     功能                      |
    | :-------------------: | :-------------------------------------------: |
    |        bin(x)         |    把十进制数x转换为二进制字符串（0b开头）    |
    |        oct(x)         |    把十进制数x转换为八进制字符串（0o开头)     |
    |        hex(x)         |   把十进制数x转换为十六进制字符串（0x开头)    |
    |       eval(str)       |               计算表达式str的值               |
    |        ord(ch)        |          返回字符ch的编码（int不行）          |
    |      round(x,y)       |   对x进行四舍五入，保留y位小数（默认y为0）    |
    |       int(st,x)       | 保留整数；把x进制下的字符串st转换为10进制的值 |
    | range(begin,end,step) |        [begin,end)公差为step的等差数列        |
    
8.  python不换行输出：

    ```python
    print(i,end='')
    ```

9.  python模块的引用：

    1.  import 模块名 as 别名（可以用别名来访问模块的对象）
    2.  from 模块名 import 对象名 as 别名

10.	python list 用法：

    |          用法          |                       作用                       |
| :--------------------: | :----------------------------------------------: |
    |     list.append(x)     |               将x添加到list的尾部                |
    |  list1.extend(list2)   |             将list2添加到list1的尾部             |
    |   list.insert(pos,x)   |            在list下标为pos的位置插入x            |
    |     list.remove(x)     |      删除list里第一次出现的x（没有会报错）       |
    |     list.pop(pos)      |          删除并返回list下标为pos的元素           |
    |     list.index(x)      |           返回list中x第一次出现的下标            |
    |     list.reverse()     |                     反转list                     |
    | list[begin​：end：​step] | 返回由list从begin到end步长为step的元素构成的列表 |
    
11.  python元组只支持通过切片操作来访问里面的元素，但如果元组的某个元素是可变的（比如列表），那么该可变元素的修改不受元组限制

12.  python 字典操作：

     ```python
     M={}
     M['name']='free'
     M['age']=18
     print(M['name'])
     #free
     print(M.get('name'))
     #free
     for i in M.items():
         print(i)
     #('name', 'free')
     #('age', 18)
     for key,value in M.items():
         print(key,value)
     #name free
     #age 18
     ```

13.  python集合操作：

     | 函数          | 作用                                                         |
     | ------------- | ------------------------------------------------------------ |
     | set.add(x)    | 在集合中插入x                                                |
     | set.pop()     | 弹出并返回集合中随机的一个的元素（有规律可循但建议不要作死地认为就是弹出首个元素） |
     | set.remove(x) | 在集合中删除x                                                |

14.  python的高级数据结构之heapq(import heapq)

     | 函数                | 作用                    |
     | ------------------- | ----------------------- |
     | heapq.heappush(H,x) | 在H中插入x              |
     | heapq.heappop(H)    | 删除并返回H中最小的元素 |

     注意：heap中的元素不是排好序的，只能保证弹出的元素是最小的

15.  python的高级数据结构之queue:

     ```python
     from queue import Queue
     q = Queue()
     q.put(1)
     q.put(2)
     q.put(3)
     print(q.qsize())
     #3
     print(q.get())
     #1
     q.put(4)
     print(q.queue)
     #deque([2, 3, 4])
     while(q.qsize()!=0):
         print(q.get())
     #2
     #3
     #4
     ```

16.  python 的高级数据结构之PriorityQueue:

     ```python
     from queue import PriorityQueue
     q = PriorityQueue()
     q.put(3)
     q.put(2)
     q.put(1)
     print(q.qsize())
     #3
     print(q.get())
     #1
     q.put(4)
     print(q.queue)
     #[2, 3, 4]
     while(q.qsize()!=0):
         print(q.get())
     #2
     #3
     #4
     ```

17.  python的高级数据结构之stack（list就可以实现）:

     ```python
     s=[1,2,3]
     s.append(4)
     print(s)
     #[1, 2, 3, 4]
     s.pop()
     print(s.pop())
     #3
     print(type(s))
     #<class 'list'>
     ```

18.  python里for和while循环后加else表示当循环正常结束（不是因为break导致结束）后运行的语句。

19.  python是基于值的内存管理方式，如下

     ```python
     for i in range(5):
         print(i)
         if i%2==0:
             i+=1
     #0
     #1
     #2
     #3
     #4
     ```

     这个样例中之所以`i+=1`没用起到任何作用是因为每次i的值发生变化时，其内存地址也会发生变化，所以每次进入循环时的i已经不再是上一次的循环变量了。

20.  python的math库：

     | 名称         | 作用                                  |
     | ------------ | ------------------------------------- |
     | sqrt()       | 求平方根                              |
     | ceil()       | 取整，和int()差不多                   |
     | gcd()        | 最大公约数                            |
     | sin()/asin() | 三角与反三角函数                      |
     | log(a,b)     | 求以a为底b的对数(也自带了log2与log10) |
     | factorial()  | 求阶乘                                |




