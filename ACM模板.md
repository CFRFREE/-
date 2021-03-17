### 2-sat

```c++
#include <bits/stdc++.h>
#define INF 2147483647
#define LL long long
#define N 2000005
using namespace std;
int n, m, cnt, dfn[N], low[N], in[N], tot;
int to[N], id[N], next[N], from[N];
stack<int> S;
inline int read()
{
    int X = 0, w = 0;
    char ch = 0;
    while (!isdigit(ch))
    {
        w |= ch == '-';
        ch = getchar();
    }
    while (isdigit(ch))
        X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
    return w ? -X : X;
}
inline void add(int x, int y)
{
    to[+cnt] = y;
    next[cnt] = from[x];
    from[x] = cnt;
}
void tarjan(int x)
{
    dfn[x] = low[x] = ++cnt;
    S.push(x);
    in[x] = 1;
    for (int i = from[x]; i; i = next[i])
    {
        int y = to[i];
        if (!dfn[y])
        {
            tarjan(y);
            low[x] = min(low[x], low[y]);
        }
        else if (in[y])
            low[x] = min(low[x], dfn[y]);
    }
    if (low[x] == dfn[x])
    {
        tot++;
        while (S.top() != x)
        {
            int X = S.top();
            id[X] = tot;
            in[X] = 0;
            S.pop();
        }
        id[S.top()] = tot;
        in[S.top()] = 0;
        S.pop();
    }
}
int main()
{
    n = read(), m = read();
    while (m--)
    {
        int a = read();
        int x = read();
        int b = read();
        int y = read();//a为x或b为y
        if (x == 0 && y == 0)
        {
            add(a + n, b);
            add(b + n, a);
        }
        if (x == 0 && y == 1)
        {
            add(a + n, b + n);
            add(b, a);
        }
        if (x == 1 && y == 0)
        {
            add(a, b);
            add(b + n, a + n);
        }
        if (x == 1 && y == 1)
        {
            add(a, b + n);
            add(b, a + n);
        }
    }
    cnt = 0;
    for (int i = 1; i <= n << 1; i++)
        if (!dfn[i])
            tarjan(i);
    for (int i = 1; i <= n; i++)
        if (id[i] == id[i + n])
        {
            printf("IMPOSSIBLE\n");
            return 0;
        }
    printf("POSSIBLE\n");
    for (int i = 1; i <= n; i++)
        if (id[i] > id[i + n])
            printf("1 ");
        else
            printf("0 ");
    return 0;
}
```

********

### kmp

```c++
#include<bits/stdc++.h>
#define LL long long
#define maxn 1000005
using namespace std;

int next[maxn], ans, len1, len2;
string st1, st2;
vector<int>pos;

inline int read()
{
	int X = 0, w = 0; char ch = 0;
	while (!isdigit(ch)) {w |= ch == '-'; ch = getchar();}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}

void KMP()
{
	next[0] = -1;
	int k = -1;
	for (int q = 1; q < len2; q++)
	{
		while (k != -1 && st2[k + 1] != st2[q])
			k = next[k];
		if (st2[k + 1] == st2[q])
			k++;
		next[q] = k;
	}

	k = -1;
	for (int i = 0; i < len1; i++)
	{
		while (k != -1 && st2[k + 1] != st1[i])
			k = next[k];
		if (st2[k + 1] == st1[i]) k++;
		if (k == len2 - 1)
		{
			i = i - len2 + 1;
			k = -1;
			cout << i + 1 << endl;//匹配成功
		}
	}
}

int main()
{
	cin >> st1;
	cin >> st2;
	len1 = st1.length();
	len2 = st2.length();
	KMP();
	for (int i = 0; i < len2; i++)
		cout << next[i] + 1 << " ";
	return 0;
}

```

******

### st表

```c++
#include<bits/stdc++.h>
using namespace std;
int n, m, st[1000001][20];
inline int read()
{
    int X = 0, w = 0; char ch = 0;
    while (!isdigit(ch)) {w |= ch == '-'; ch = getchar();}
    while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
    return w ? -X : X;
}
inline int max(int x, int y) {
    return x > y ? x : y;
}
inline int search(int x, int y)
{
    int k = log(y - x + 1) / log(2);
    return max(st[x][k], st[y - (1 << k) + 1][k]);
}
int main()
{
    n = read(), m = read();
    for (register int i = 1; i <= n; i++)
        st[i][0] = read();
    for (register int j = 1; j <= 20; j++)
        for (register int i = 1; i + (1 << j) - 1 <= n; i++)
            st[i][j] = max(st[i][j - 1], st[i + (1 << j - 1)][j - 1]);
    for (register int i = 1; i <= m; i++)
    {
        int x = read(), y = read();
        printf("%d\n", search(x, y));
    }
    return 0;
}
```

*****

### tarjan LCA

```c++
#include<bits/stdc++.h>
#define LL long long
#define maxn 500005

using namespace std;

int n, m, s, vis[maxn], ans[maxn], fa[maxn];
vector<int>v[maxn];
vector< pair<int, int> >ask[maxn];

inline int read()
{
    int X = 0, w = 0; char ch = 0;
    while (!isdigit(ch)) {w |= ch == '-'; ch = getchar();}
    while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
    return w ? -X : X;
}

int find(int x)
{
    if (fa[x] != x)fa[x] = find(fa[x]);
    return fa[x];
}

void tarjan(int x)
{
    vis[x] = 1;
    int X = find(x);
    for (register int i = 0; i < v[x].size(); i++)
    {
        int y = v[x][i];
        if (vis[y])continue;
        tarjan(y);
        int Y = find(y);
        if (X != Y)fa[Y] = X;
    }

    for (register int i = 0; i < ask[x].size(); i++)
    {
        int y = ask[x][i].first, id = ask[x][i].second;
        if (vis[y])ans[id] = find(y);
    }
}

int main()
{
    n = read(), m = read(), s = read();

    for (register int i = 1; i < n; i++)
    {
        int x = read(), y = read();
        v[x].push_back(y);
        v[y].push_back(x);
    }

    for (register int i = 1; i <= m; i++)
    {
        int x = read(), y = read();
        ask[x].push_back(make_pair(y, i));
        ask[y].push_back(make_pair(x, i));
    }

    for (register int i = 1; i <= n; i++)
        fa[i] = i;

    tarjan(s);

    for (register int i = 1; i <= m; i++)
        printf("%d\n", ans[i]);

    return 0;
}
```

*****

### tarjan缩点

```c++
#include<bits/stdc++.h>
#define ll long long
using namespace std;
const int maxn = 100005;
vector<int>a[maxn];
stack<int>s;
queue<int>q;
int n, m, vis[maxn], v[maxn], value[maxn], from[maxn], to[maxn], sum;
int ans, dfn[maxn], low[maxn], cnt, in[maxn], id[maxn], d[maxn], N, f[maxn];
inline int read()
{
	int X = 0, w = 0; char ch = 0;
	while (!isdigit(ch)) {w |= ch == '-'; ch = getchar();}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}
inline void tarjan(int x)
{
	dfn[x] = low[x] = ++cnt;
	s.push(x);
	in[x] = 1;
	if (a[x].size())
		for (int i = 0; i < a[x].size(); i++)
		{
			if (!dfn[ a[x][i] ])
			{
				tarjan(a[x][i]);
				low[x] = min(low[x], low[ a[x][i] ]);
			}
			else if (in[ a[x][i] ]) {
				low[x] = min(low[x], dfn[ a[x][i] ]);
			}
		}
	if (dfn[x] == low[x])
	{
		N++;
		while (s.top() != x)
		{
			int X = s.top();
			s.pop();
			id[X] = N;
			in[X] = 0;
			value[N] += v[X];
		}
		in[s.top()] = 0;
		id[s.top()] = N;
		value[N] += v[s.top()];
		s.pop();
	}
}
inline void search(int x)
{
	if (f[x])return;
	f[x] = value[x];
	int maxsum = 0;
	for (int i = 0; i < a[x].size(); i++)
	{
		if (!f[ a[x][i] ])
			search(a[x][i]);
		maxsum = max(maxsum, f[ a[x][i] ]);
	}
	f[x] += maxsum;
}
int main()
{
	n = read(), m = read();
	for (int i = 1; i <= n; i++)
		v[i] = read();
	for (int i = 1; i <= m; i++)
	{
		from[i] = read(), to[i] = read();
		a[from[i]].push_back(to[i]);
	}
	for (int i = 1; i <= n; i++)
	{
		if (!dfn[i])tarjan(i);
	}
	for (int i = 0; i < maxn; i++)
		a[i].clear();
	for (int i = 1; i <= m; i++)
		if (id[ from[i] ] != id[ to[i] ])
			a[ id[ from[i] ] ].push_back(id[ to[i] ]);
	for (int i = 1; i <= N; i++)
		search(i);//spfa(i);
	for (int i = 1; i <= N; i++)
		ans = max(ans, f[i]);
	cout << ans;
	return 0;
}
```

****

### 倍增LCA

```c++
#include<bits/stdc++.h>
#define LL long long
using namespace std;
const int maxn = 500001;
int n, m, s, lg[maxn], deep[maxn], fa[maxn][22];
vector<int>v[maxn];

inline int read()
{
	int X = 0, w = 0; char ch = 0;
	while (!isdigit(ch)) {w |= ch == '-'; ch = getchar();}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}

inline void dfs(int x, int father)
{
	deep[x] = deep[father] + 1;
	fa[x][0] = father;
	for (int i = 1; (1 << i) <= deep[x]; i++)
		fa[x][i] = fa[fa[x][i - 1]][i - 1];
	for (int i = 0; i < v[x].size(); i++)
		if (v[x][i] != father)dfs(v[x][i], x);
}

inline int LCA(int x, int y)
{
	if (deep[x] < deep[y])swap(x, y);
	while (deep[x] != deep[y])x = fa[x][lg[deep[x] - deep[y]] - 1];
	if (x == y)return x;
	for (register int i = lg[deep[x]]; i >= 0; i--)
		if (fa[x][i] != fa[y][i])x = fa[x][i], y = fa[y][i];
	return fa[x][0];
}

int main()
{
	n = read(), m = read(), s = read();
	for (register int i = 1; i <= n; i++)
	{
		lg[i] = lg[i - 1];
		if (i == (1 << lg[i - 1]))lg[i]++;
	}
	for (register int i = 1; i <= n - 1; i++)
	{
		int x = read(), y = read();
		v[x].push_back(y);
		v[y].push_back(x);
	}

	dfs(s, 0);
	for (register int i = 1; i <= m; i++)
	{
		int x = read(), y = read();
		printf("%d\n", LCA(x, y));
	}

	return 0;
}
```

****

### 并查集

```c++
#include <bits/stdc++.h>
#define INF 2147483647
#define LL long long
#define N 1000005

using namespace std;
int n, m, fa[N];
inline int read()
{
	int X = 0, w = 0;
	char ch = 0;
	while (!isdigit(ch))
	{
		w |= ch == '-';
		ch = getchar();
	}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}
int find(int x)
{
	if (x == fa[x])return x;
	fa[x] = find(fa[x]);
	return fa[x];
}
int main()
{
	n = read(), m = read();
	for (int i = 1; i <= n; i++)
		fa[i] = i;
	for (int i = 1; i <= m; i++)
	{
		int z = read(), x = read(), y = read();
		if (z == 1)
			fa[find(x)] = fa[find(y)];
		else if (find(fa[x]) == find(fa[y]))
			printf("Y\n");
		else printf("N\n");
	}
	return 0;
}
```

*****

### 差分约束

```c++
#include <bits/stdc++.h>
#define INF 2147483647
#define LL long long
#define N 100005
using namespace std;
int n, m, cnt, sum[N], from[N], to[N], next[N], vis[N], d[N], v[N];
queue<int>q;
inline int read()
{
	int X = 0, w = 0;
	char ch = 0;
	while (!isdigit(ch))
	{
		w |= ch == '-';
		ch = getchar();
	}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}
inline void add(int x, int y, int z)
{
	to[++cnt] = y;
	v[cnt] = z;
	next[cnt] = from[x];
	from[x] = cnt;
}
int main()
{
	n = read(), m = read();
	for (int i = 1; i <= m; i++)
	{
		int x = read(), y = read(), z = read();
		add(y, x, z);
	}
	for (int i = 1; i <= n; i++)
		add(0, i, 0);
	memset(d, 0x3f, sizeof(d));
	d[0] = 0;
	vis[0] = 1;
	q.push(0);
	while (!q.empty())
	{
		int x = q.front();
		q.pop();
		vis[x] = 0;
		for (int i = from[x]; i; i = next[i])
		{
			int y = to[i];
			if (d[y] > d[x] + v[i])
			{
				d[y] = d[x] + v[i];
				if (!vis[y])
				{
					sum[y]++;
					q.push(y);
					if (sum[y] > n)
					{
						printf("NO\n");
						return 0;
					}
					vis[y] = 1;
				}
			}
		}
	}
	for (int i = 1; i <= n; ++i)
		printf("%d ", d[i]);
	return 0;
}
```

*****

### 堆优化dij

````cpp
#include<bits/stdc++.h>
#define LL long long
#define MAXN 100010

using namespace std;
int n, m, s, cnt, d[MAXN], vis[MAXN], to[MAXN * 2], from[MAXN * 2], next[MAXN * 2], v[MAXN * 2];

void add(int x, int y, int z)
{
	to[++cnt] = y;
	v[cnt] = z;
	next[cnt] = from[x];
	from[x] = cnt;
}

int main()
{
	scanf("%d%d%d", &n, &m, &s);
	for (int i = 1; i <= m; i++)
	{
		int x, y, z;
		scanf("%d%d%d", &x, &y, &z);
		add(x, y, z);
	}
	memset(d, 0x3f, sizeof(d));
	d[s] = 0;
	priority_queue<pair<int, int> >q;
	q.push(make_pair(0, s));
	while (!q.empty())
	{
		int val = q.top().first, x = q.top().second;
		q.pop();
		if (val > d[x])continue;
		if (vis[x])continue;
		vis[x] = 1;
		for (int i = from[x]; i; i = next[i])
		{
			int y = to[i];
			if (d[x] + v[i] < d[y])
			{
				d[y] = d[x] + v[i];
				q.push(make_pair(-d[y], y));
			}
		}
	}
	for (int i = 1; i <= n; i++)
		printf("%d ", d[i]);
	return 0;
}
````

*****

### 二分

```cpp
#include <bits/stdc++.h>
#define INF 2147483647
#define LL long long
#define N 100005
using namespace std;
int a[N];
inline int read()
{
	int X = 0, w = 0;
	char ch = 0;
	while (!isdigit(ch))
	{
		w |= ch == '-';
		ch = getchar();
	}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}
int main()
{
	/*
	while(l<=r)
	{
		LL mid=(l+r)>>1;
		if(ok(mid))
		{
			ans=min(ans,mid);
			r=mid-1;
		}
		else l=mid+1;
	}
	*/
	int n = read(), k = read();
	for (int i = 1; i <= n; i++)
		a[i] = read();
	int pos3 = lower_bound(a + 1, a + n + 1, k, greater<int>()) - a - 1;
	printf("%d\n", pos3);
	return 0;
}
```

*****

### 二进制穷举

```c++
#include<stdio.h>
#include<string.h>
int main()
{
	int n;
	scanf("%d",&n);
	int i;
	char s[100],l;
	for(i=0;i<(1<<n);i++)
	{
		int tmp=i;
		l=0;
		while(tmp)
		{
			if(tmp&1) s[l++]='1';
			else s[l++]='0';
			tmp>>=1;
		}
		s[l]='\0';
		strrev(s);//倒序
		puts(s);
	}
	return 0;
}
```

*****

### 分块

```c++
#include<bits/stdc++.h>
#define LL long long
#define maxn 400005
using namespace std;

LL a[maxn], b[maxn], tot[maxn], add[maxn], sum[maxn], Left[maxn], Right[maxn];

inline LL read()
{
    LL X = 0, w = 0;
    char ch = 0;
    while (!isdigit(ch))
    {
        w |= ch == '-';
        ch = getchar();
    }
    while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
    return w ? -X : X;
}

int main()
{
    LL n = read(), m = read();
    LL len = sqrt(n);
    LL tot = n / len;
    if (n % tot) tot++;

    for (LL i = 1; i <= n; i++)
    {
        a[i] = read();
        b[i] = (i - 1) / len + 1;
        sum[ b[i] ] += a[i];
    }

    for (LL i = 1; i <= tot; i++)
    {
        Left[i] = (i - 1) * len + 1;
        Right[i] = i * len;
    }

    while ( m-- )
    {
        int p = read();
        if (p == 1)
        {
            LL x = read(), y = read(), z = read();

            if (b[x] == b[y])
            {
                for (int i = x; i <= y; i++)
                    a[i] += z;
                sum[ b[x] ] += z * (y - x);
            }
            else
            {
                for (LL i = b[x] + 1; i <= b[y] - 1; i++)
                    add[i] += z;

                for (LL i = x; i <= Right[ b[x] ]; i++)
                    a[i] += z, sum[ b[x] ] += z;

                for (LL i = Left[ b[y] ]; i <= y; i++)
                    a[i] += z, sum[ b[y] ] += z;
            }
        }
        else
        {
            LL ans = 0;
            LL x = read(), y = read();

            if (b[x] == b[y])
            {
                for (int i = x; i <= y; i++)
                    ans += a[i] + add[ b[i] ];
            }
            else
            {
                for (LL i = b[x] + 1; i <= b[y] - 1; i++)
                    ans += (sum[i] + (add[i] * len));

                for (LL i = x; i <= min(Right[ b[x] ], y); i++)
                    ans += a[i] + add[ b[i] ];

                for (LL i = max(Left[ b[y] ], x); i <= y; i++)
                    ans += a[i] + add[ b[i] ];
            }
            cout << ans << endl;
        }
    }
    return 0;
}
```

*****

### 割点

````c++
//删除后强连通分量数量增加
#include<bits/stdc++.h>
using namespace std;
const int maxn = 100001;
int n, m, cnt, dfn[maxn], low[maxn], father[maxn];
vector<int>v[maxn];
set<int>s;
inline int read()
{
	int X = 0, w = 0; char ch = 0;
	while (!isdigit(ch)) {w |= ch == '-'; ch = getchar();}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}
void tarjan(int x)
{
	dfn[x] = low[x] = ++cnt;
	int child = 0;
	for (int i = 0; i < v[x].size(); i++)
	{
		if (!dfn[v[x][i]])
		{
			father[v[x][i]] = x;
			child++;
			tarjan(v[x][i]);
			low[x] = min(low[x], low[v[x][i]]);

			if (father[x] && low[v[x][i]] >= dfn[x])
				s.insert(x);
		} else if (father[v[x][i]] != x)
			low[x] = min(low[x], dfn[v[x][i]]);
	}
	if (!father[x] && child >= 2)
		s.insert(x);
}
int main()
{
	n = read(), m = read();
	for (int i = 1; i <= m; i++)
	{
		int x = read(), y = read();
		v[x].push_back(y);
		v[y].push_back(x);
	}

	for (int i = 1; i <= n; i++)
		if (!dfn[i])
		{
			//	cnt=0;
			tarjan(i);
		}

	cout << s.size() << endl;
	set<int>::iterator it;
	for (it = s.begin(); it != s.end(); it++)
		cout << *it << " ";
	return 0;
}
````

*****

### 莫队

```c++
//袜子颜色的概率
#include<bits/stdc++.h>
#define LL long long
#define MAXN 100005

using namespace std;

LL n,m,len,tot,a[MAXN],b[MAXN],cnt[MAXN];
LL answer,ans_fz[MAXN],ans_fm[MAXN];

struct Node
{
    LL l,r,id;
} ask[MAXN];

inline LL read()
{
    LL X=0,w=0;
    char ch=0;
    while(!isdigit(ch))
    {
        w|=ch=='-';
        ch=getchar();
    }
    while(isdigit(ch)) X=(X<<3)+(X<<1)+(ch^48),ch=getchar();
    return w?-X:X;
}

inline LL cmp(Node x,Node y)
{
    if(b[x.l]!=b[y.l])return x.l<y.l;
    if(b[x.l]%2==1)return x.r<y.r;
    return x.r>y.r;
}

inline void add(LL pos)
{
    answer+=cnt[a[pos]]*2+1;
    ++cnt[a[pos]];
}

inline void remove(LL pos)
{
    answer+=1-cnt[a[pos]]*2;
    --cnt[a[pos]];
}

int main()
{
    n=read(),m=read();
    len=int(pow(n,0.6));
    tot=n/len;
    if(n%len)tot++;
    for(register LL i=1; i<=n; i++)
    {
        a[i]=read();
        b[i]=(i-1)/len+1;
    }
    for(register LL i=1; i<=m; i++)
    {
        ask[i].l=read();
        ask[i].r=read();
        ask[i].id=i;
        ans_fm[i]=1;
    }

    sort(ask+1,ask+m+1,cmp);

    LL curl=ask[1].l,curr=ask[1].r;

    for(register LL i=curl; i<=curr; i++)
        add(i);

    if(answer!=ask[1].r-ask[1].l+1)
    {
        LL fz=answer-ask[1].r+ask[1].l-1;
        LL fm=(ask[1].r-ask[1].l+1)*(ask[1].r-ask[1].l);
        LL g=__gcd(fz,fm);
        ans_fz[ask[1].id]=fz/g;
        ans_fm[ask[1].id]=fm/g;
    }

    for(register LL i=2; i<=m; i++)
    {
        while(curl<ask[i].l)remove(curl++);
        while(curl>ask[i].l)add(--curl);
        while(curr>ask[i].r)remove(curr--);
        while(curr<ask[i].r)add(++curr);
        if(answer!=ask[i].r-ask[i].l+1)
        {
            LL fz=answer-ask[i].r+ask[i].l-1;
            LL fm=(ask[i].r-ask[i].l+1)*(ask[i].r-ask[i].l);
            LL g=__gcd(fz,fm);
            ans_fz[ask[i].id]=fz/g;
            ans_fm[ask[i].id]=fm/g;
        }
    }

    for(register LL i=1; i<=m; i++)
        printf("%d/%d\n",ans_fz[i],ans_fm[i]);

    return 0;
}
```

******

### Fhq Treap

```c++
// luogu-judger-enable-o2
#include<bits/stdc++.h>
#define LL long long
#define MAXN 100005

using namespace std;
int n, cnt, root;
struct Node {
	int v, rnd, son[2], size;
} tree[MAXN];

inline int read()
{
	int X = 0, w = 0;
	char ch = 0;
	while (!isdigit(ch))
	{
		w |= ch == '-';
		ch = getchar();
	}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}

void PushUp(int k)
{
	tree[k].size = tree[tree[k].son[0]].size + tree[tree[k].son[1]].size + 1;
}

int build(int v)
{
	++cnt;
	tree[cnt].size = 1;
	tree[cnt].v = v;
	tree[cnt].rnd = rand();
	return cnt;
}

int merge(int x, int y)
{
	if (!x || !y)return x + y;
	if (tree[x].rnd < tree[y].rnd)
	{
		tree[x].son[1] = merge(tree[x].son[1], y);
		PushUp(x);
		return x;
	}
	else
	{
		tree[y].son[0] = merge(x, tree[y].son[0]);
		PushUp(y);
		return y;
	}
}

void split(int now, int k, int &x, int &y)
{
	if (!now)x = y = 0;
	else
	{
		if (tree[now].v <= k)
		{
			x = now;
			split(tree[now].son[1], k, tree[now].son[1], y);
		}
		else
		{
			y = now;
			split(tree[now].son[0], k, x, tree[now].son[0]);
		}
		PushUp(now);
	}
}

int kth(int now, int k)
{
	while (1)
	{
		if (k <= tree[tree[now].son[0]].size)
			now = tree[now].son[0];
		else
		{
			if (k == tree[tree[now].son[0]].size + 1)
				return now;
			else
			{
				k -= tree[tree[now].son[0]].size + 1;
				now = tree[now].son[1];
			}
		}
	}
}

int main()
{
	srand(time(NULL));
	n = read();
	int x, y, z;
	while (n--)
	{
		int opt = read(), a = read();
		if (opt == 1)
		{
			split(root, a, x, y);
			root = merge(merge(x, build(a)), y);
		}//插入x
		if (opt == 2)
		{
			split(root, a, x, z);
			split(x, a - 1, x, y);
			y = merge(tree[y].son[0], tree[y].son[1]);
			root = merge(merge(x, y), z);
		}//删除x
		if (opt == 3)
		{
			split(root, a - 1, x, y);
			printf("%d\n", tree[x].size + 1);
			root = merge(x, y);
		}//查询x的排名
		if (opt == 4)
			printf("%d\n", tree[kth(root, a)].v);//查询排名为x的数
		if (opt == 5)
		{
			split(root, a - 1, x, y);
			printf("%d\n", tree[kth(x, tree[x].size)].v);
			root = merge(x, y);
		}//求x前驱
		if (opt == 6)
		{
			split(root, a, x, y);
			printf("%d\n", tree[kth(y, 1)].v);
			root = merge(x, y);
		}//求x后继
	}
	return 0;
}
```

*****

### 全排列

```cpp
#include<bits/stdc++.h>
#define MAXN 50005

using namespace std;

int n, m, d[6][MAXN], vis[MAXN], last[200005], next[200005], to[200005], value[200005],
    home[6], sum, ans = 1e18, cnt;
queue<int>q;
map<int, int>M;

inline int read()
{
	int X = 0, w = 0;
	char ch = 0;
	while (!isdigit(ch))
	{
		w |= ch == '-';
		ch = getchar();
	}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}

inline void add(int u, int v, int w)
{
	to[++cnt] = v;
	next[cnt] = last[u];
	last[u] = cnt;
	value[cnt] = w;
}

void spfa(int from, int id)
{
	memset(d[id], 0x7f, sizeof(d[id]));
	memset(vis, 0, sizeof(vis));
	d[id][from] = 0;
	vis[from] = 1;
	q.push(from);
	while (!q.empty())
	{
		int x = q.front();
		q.pop();
		vis[x] = 0;
		for (register int i = last[x]; i; i = next[i])
		{
			int y = to[i];
			if (d[id][x] + value[i] < d[id][y])
			{
				d[id][y] = d[id][x] + value[i];
				if (!vis[y])
				{
					vis[y] = 1;
					q.push(y);
				}
			}
		}
	}
}

int main()
{
	n = read(), m = read();
	for (register int i = 1; i <= 5; i++)
		home[i] = read();

	while (m--)
	{
		int x = read(), y = read(), z = read();
		add(x, y, z);
		add(y, x, z);
	}

	sort(home + 1, home + 6);
	for (register int i = 1; i <= 5; i++)
	{
		spfa(home[i], i);
		M[home[i]] = i;
	}

	do
	{
		int sum = d[M[home[1]]][1];
		for (register int i = 1; i <= 4; i++)
			sum += d[M[home[i]]][home[i + 1]];
		if (sum < ans)ans = sum;
	}
	while (next_permutation(home + 1, home + 6));

	printf("%d", ans);
	return 0;
}
```

******

### 快读

```cpp
#include<bits/stdc++.h>
#define LL long long
#define N 100005
using namespace std;
inline int read()
{
	int X = 0, w = 0;
	char ch = 0;
	while (!isdigit(ch))
	{
		w |= ch == '-';
		ch = getchar();
	}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}
int main()
{

	return 0;
}
```

*****

*****

### 树链剖分

```cpp
// luogu-judger-enable-o2
#include<bits/stdc++.h>
#define LL long long
#define maxn 100005
using namespace std;
int n, m, r, P, id[maxn], fa[maxn], son[maxn], cnt;
int top[maxn], sum, a[maxn], b[maxn], deep[maxn], size[maxn];
vector<int>v[maxn];
struct node {
	int l, r, w, f;
} tree[maxn << 2];
inline int read()
{
	int X = 0, w = 0; char ch = 0;
	while (!isdigit(ch)) {w |= ch == '-'; ch = getchar();}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}

void PushUp(int k)
{
	tree[k].w = tree[k << 1].w + tree[k << 1 | 1].w;
	tree[k].w %= P;
}

void PushDown(int k)
{
	if (tree[k].f)
	{
		int x = tree[k].f;
		tree[k].f = 0;
		tree[k << 1].f += x;
		tree[k << 1 | 1].f += x;
		tree[k << 1].w += x * (tree[k << 1].r - tree[k << 1].l + 1);
		tree[k << 1 | 1].w += x * (tree[k << 1 | 1].r - tree[k << 1 | 1].l + 1);
		tree[k << 1].w %= P;
		tree[k << 1 | 1].w %= P;
	}
}

void build(int k, int ll, int rr)
{
	tree[k].l = ll;
	tree[k].r = rr;
	if (ll == rr)
	{
		tree[k].w = a[ll];
		return;
	}
	int mid = (ll + rr) >> 1;
	build(k << 1, ll, mid);
	build(k << 1 | 1, mid + 1, rr);
	PushUp(k);
}

void change(int k, int ll, int rr, int d)
{
	if (tree[k].l >= ll && tree[k].r <= rr)
	{
		tree[k].f += d;
		tree[k].w += d * (tree[k].r - tree[k].l + 1);
		return;
	}
	PushDown(k);
	int mid = (tree[k].l + tree[k].r) >> 1;
	if (mid >= ll)change(k << 1, ll, rr, d);
	if (rr > mid)change(k << 1 | 1, ll, rr, d);
	PushUp(k);
}

int search(int k, int ll, int rr)
{
	int ans = 0;
	if (tree[k].l >= ll && tree[k].r <= rr)
		return tree[k].w % P;

	PushDown(k);
	int mid = (tree[k].l + tree[k].r) >> 1;
	if (mid >= ll)ans += search(k << 1, ll, rr);
	if (rr > mid)ans += search(k << 1 | 1, ll, rr);
	return ans % P;
}

int dfs1(int x, int f, int dep)
{
	deep[x] = dep;
	fa[x] = f;
	size[x] = 1;
	int maxson = -1;
	for (int i = 0; i < v[x].size(); i++)
	{
		int y = v[x][i];
		if (y == fa[x])continue;
		size[x] += dfs1(y, x, dep + 1);
		if (size[y] > maxson)
		{
			son[x] = y;
			maxson = size[y];
		}
	}
	return size[x];
}

void dfs2(int x, int topf)
{
	top[x] = topf;
	id[x] = ++cnt;
	a[cnt] = b[x];

	if (!son[x])return;

	dfs2(son[x], topf);

	for (int i = 0; i < v[x].size(); i++)
		if (!id[v[x][i]])
			dfs2(v[x][i], v[x][i]);
}

void path_change(int x, int y, int z)
{
	while (top[x] != top[y])
	{
		if (deep[top[x]] < deep[top[y]])swap(x, y);
		change(1, id[top[x]], id[x], z);
		x = fa[top[x]];
	}
	if (deep[x] > deep[y])swap(x, y);
	change(1, id[x], id[y], z);
}

void path_search(int x, int y)
{
	int ans = 0;
	while (top[x] != top[y])
	{
		if (deep[top[x]] < deep[top[y]])swap(x, y);
		ans += search(1, id[top[x]], id[x]);
		ans %= P;
		x = fa[top[x]];
	}
	if (deep[x] > deep[y])swap(x, y);
	ans += search(1, id[x], id[y]);
	cout << ans % P << endl;
}

void tree_change(int x, int y)
{
	change(1, id[x], id[x] + size[x] - 1, y);
}

void tree_search(int x)
{
	cout << search(1, id[x], id[x] + size[x] - 1) << endl;
}

int main()
{
	n = read(), m = read(), r = read(), P = read();
	for (int i = 1; i <= n; i++)
		b[i] = read();
	for (int i = 1; i < n; i++)
	{
		int x = read(), y = read();
		v[x].push_back(y);
		v[y].push_back(x);
	}

	dfs1(r, 0, 1);

	dfs2(r, r);
	build(1, 1, n);

	for (int i = 1; i <= m; i++)
	{
		int p = read();
		if (p == 1)
		{
			int x = read(), y = read(), z = read();
			path_change(x, y, z);
		}
		if (p == 2)
		{
			int x = read(), y = read();
			path_search(x, y);
		}
		if (p == 3)
		{
			int x = read(), y = read();
			tree_change(x, y);
		}
		if (p == 4)
		{
			int x = read();
			tree_search(x);
		}
	}
	return 0;
}
```

*****

### 树状数组

```cpp
#include <bits/stdc++.h>
#define INF 2147483647
#define LL long long
#define N 500005

using namespace std;
int n, m, tree[N];
inline int read()
{
	int X = 0, w = 0;
	char ch = 0;
	while (!isdigit(ch))
	{
		w |= ch == '-';
		ch = getchar();
	}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}
inline int low(int x)
{
	return x & -x;
}
void add(int x, int k)
{
	while (x <= n)
	{
		tree[x] += k;
		x += low(x);
	}
}
LL sum(int x)
{
	LL s = 0;
	while (x != 0)
	{
		s += tree[x];
		x -= low(x);
	}
	return s;
}
int main()
{
	n = read(), m = read();
	for (int i = 1; i <= n; i++)
	{
		int X = read();
		add(i, X);
	}
	for (int i = 0; i < m; ++i)
	{
		int p = read(), t = read(), k = read();
		if (p == 1)
			add(t, k);
		else
			printf("%lld\n", sum(k) - sum(t - 1));
	}
	return 0;
}
```

### 线段树

```c++
#include <bits/stdc++.h>
#define INF 2147483647
#define LL long long
#define N 1000005

using namespace std;
int n, m;
struct Node {
	LL l, r, f, w;
} tree[N << 2];
inline int read()
{
	int X = 0, w = 0;
	char ch = 0;
	while (!isdigit(ch))
	{
		w |= ch == '-';
		ch = getchar();
	}
	while (isdigit(ch)) X = (X << 3) + (X << 1) + (ch ^ 48), ch = getchar();
	return w ? -X : X;
}
inline void PushUp(int k)
{
	tree[k].w = tree[k << 1].w + tree[k << 1 | 1].w;
}
inline void PushDown(int k)
{
	if (!tree[k].f)return;
	LL x = tree[k].f;
	tree[k].f = 0;
	tree[k << 1].f += x;
	tree[k << 1 | 1].f += x;
	tree[k << 1].w += (tree[k << 1].r - tree[k << 1].l + 1) * x;
	tree[k << 1 | 1].w += (tree[k << 1 | 1].r - tree[k << 1 | 1].l + 1) * x;

}
void build(int k, int l, int r)
{
	tree[k].l = l;
	tree[k].r = r;
	if (l == r)
	{
		tree[k].w = read();
		return;
	}
	int mid = (l + r) >> 1;
	build(k << 1, l, mid);
	build(k << 1 | 1, mid + 1, r);
	PushUp(k);
}
void change(int k, int l, int r, int d)
{
	if (tree[k].l >= l && tree[k].r <= r)
	{
		tree[k].f += d;
		tree[k].w += d * (tree[k].r - tree[k].l + 1);
		return;
	}
	PushDown(k);
	int mid = (tree[k].l + tree[k].r) >> 1;
	if (mid >= l)change(k << 1, l, r, d);
	if (mid < r)change(k << 1 | 1, l, r, d);
	PushUp(k);
}

LL query(int k, int l, int r)
{
	if (tree[k].l >= l && tree[k].r <= r)
		return tree[k].w;
	PushDown(k);
	int mid = (tree[k].l + tree[k].r) >> 1;
	LL sum = 0;
	if (mid >= l)sum += query(k << 1, l, r);
	if (mid < r)sum += query(k << 1 | 1, l, r);
	return sum;
}
int main()
{
	n = read(), m = read();
	build(1, 1, n);
	while (m--)
	{
		int opt = read();
		if (opt == 1)
		{
			int x = read(), y = read(), z = read();
			change(1, x, y, z);
		}
		else
		{
			int x = read(), y = read();
			printf("%lld\n", query(1, x, y));
		}
	}
	return 0;
}
```

******

### 匈牙利

```cpp
// luogu-judger-enable-o2
#include<bits/stdc++.h>
using namespace std;
const int maxn = 1005;
int n, m, e, match[maxn], vis[maxn], x, y, ans;
vector<int>v[maxn];
int find(int x)
{
    for (int i = 0; i < v[x].size(); i++)
    {
        int y = v[x][i];
        if (!vis[y])
        {
            vis[y] = 1;
            if (!match[y] || find(match[y])) {
                match[y] = x;
                return 1;
            }
        }
    }
    return 0;
}
int main()
{
    scanf("%d%d%d", &n, &m, &e);
    for (int i = 1; i <= e; i++)
    {
        scanf("%d%d", &x, &y);
        if (x <= n && y <= m)
            v[x].push_back(y);
    }
    for (int i = 1; i <= n; i++)
    {
        memset(vis, 0, sizeof(vis));
        if (find(i))ans++;
    }
    cout << ans;
}
```

*****

### 主席树

```cpp
//静态区间第k小
#include<bits/stdc++.h>
#define LL long long
#define maxn 200005

using namespace std;

int n, m, q, cnt;
int L[maxn << 5], R[maxn << 5], sum[maxn << 5], T[maxn], a[maxn], b[maxn];

inline int build(int l, int r)
{
    int rt = ++cnt;
    sum[rt] = 0;
    if (l < r)
    {
        int mid = (l + r) >> 1;
        L[rt] = build(l, mid);
        R[rt] = build(mid + 1, r);
    }
    return rt;
}

inline int change(int pre, int l, int r, int x)
{
    int rt = ++cnt;
    L[rt] = L[pre], R[rt] = R[pre], sum[rt] = sum[pre] + 1;
    if (l < r)
    {
        int mid = (l + r) >> 1;
        if (x <= mid)L[rt] = change(L[pre], l, mid, x);
        else R[rt] = change(R[pre], mid + 1, r, x);
    }
    return rt;
}

inline int search(int u, int v, int l, int r, int k)
{
    if (l >= r)return l;
    int x = sum[L[v]] - sum[L[u]];
    int mid = (l + r) >> 1;
    if (x >= k)return search(L[u], L[v], l, mid, k);
    else return search(R[u], R[v], mid + 1, r, k - x);
}

int main()
{
    scanf("%d%d", &n, &q);
    for (int i = 1; i <= n; i++)
    {
        scanf("%d", &a[i]);
        b[i] = a[i];
    }
    sort(b + 1, b + n + 1);
    m = unique(b + 1, b + n + 1) - b - 1;
    T[0] = build(1, m);

    for (int i = 1; i <= n; i++)
    {
        int t = lower_bound(b + 1, b + m + 1, a[i]) - b;
        T[i] = change(T[i - 1], 1, m, t);
    }

    while (q--)
    {
        int x, y, z;
        scanf("%d%d%d", &x, &y, &z);
        int t = search(T[x - 1], T[y], 1, m, z);
        printf("%d\n", b[t]);
    }

    return 0;
}
```

*****

### 字符串

```c++
string s1 = "this is ok";
string s2 = s1.substr(2, 4);  // s2 = "is i"
s2 = s1.substr(2);  // s2 = "is is ok"

#include <iostream>
#include <string>
using namespace std;
int main()
{
    string s1("Source Code");
    int n;
    if ((n = s1.find('u')) != string::npos) //查找 u 出现的位置
        cout << "1) " << n << "," << s1.substr(n) << endl;
    //输出 l)2,urce Code
    if ((n = s1.find("Source", 3)) == string::npos)
        //从下标3开始查找"Source"，找不到
        cout << "2) " << "Not Found" << endl;  //输出 2) Not Found
    if ((n = s1.find("Co")) != string::npos)
        //查找子串"Co"。能找到，返回"Co"的位置
        cout << "3) " << n << ", " << s1.substr(n) << endl;
    //输出 3) 7, Code
    if ((n = s1.find_first_of("ceo")) != string::npos)
        //查找第一次出现或 'c'、'e'或'o'的位置
        cout << "4) " << n << ", " << s1.substr(n) << endl;
    //输出 4) l, ource Code
    if ((n = s1.find_last_of('e')) != string::npos)
        //查找最后一个 'e' 的位置
        cout << "5) " << n << ", " << s1.substr(n) << endl;  //输出 5) 10, e
    if ((n = s1.find_first_not_of("eou", 1)) != string::npos)
        //从下标1开始查找第一次出现非 'e'、'o' 或 'u' 字符的位置
        cout << "6) " << n << ", " << s1.substr(n) << endl;
    //输出 6) 3, rce Code
    return 0;
}



string s1("Real Steel");
s1.replace(1, 3, "123456", 2, 4);  //用 "123456" 的子串(2,4) 替换 s1 的子串(1,3)
cout << s1 << endl;  //输出 R3456 Steel
string s2("Harry Potter");
s2.replace(2, 3, 5, '0');  //用 5 个 '0' 替换子串(2,3)
cout << s2 << endl;  //输出 HaOOOOO Potter
int n = s2.find("OOOOO");  //查找子串 "00000" 的位置，n=2
s2.replace(n, 5, "XXX");  //将子串(n,5)替换为"XXX"
cout << s2 < < endl;  //输出 HaXXX Potter


string s1("Real Steel");
s1.erase(1, 3);  //删除子串(1, 3)，此后 s1 = "R Steel"
s1.erase(5);  //删除下标5及其后面的所有字符，此后 s1 = "R Ste"

string s1("Limitless"), s2("00");
s1.insert(2, "123");  //在下标 2 处插入字符串"123"，s1 = "Li123mitless"
s1.insert(3, s2);  //在下标 2 处插入 s2 , s1 = "Li10023mitless"
s1.insert(3, 5, 'X');  //在下标 3 处插入 5 个 'X'，s1 = "Li1XXXXX0023mitless"
```