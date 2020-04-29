## G-thinker+
G-thinker+ is an improved version of G-thinker https://yanlab19870714.github.io/yanda/gthinker/

We made the following changes to accommodate "big" tasks that needs prioritized concurrent processing.

Firstly, we added a global queue shared by all computing threads to accommodate "big" tasks, so that they can be fetched by available computing threads ASAP.

Secondly, work stealing now only balances big tasks.

Thirdly, a new user defined function is added to decide whether a task is a "big" task:

bool Comper::is_bigtask(TaskT * task)

Finally, we added a new application under folder "app_quasiclique" for mining maxial quasicliques (a postprocessing is needed to remove non-maximal ones).

### What is special about G-thinker and G-thinker+?
Existing Big Data frameworks are dominantly **data-intensive** where communication caused by data movement is the bottleneck, and CPU cores are underutilized. Solving **compute-intensive** graph mining problems using these frameworks often results in even worse performance.

G-thinker fundamentally changes the data-intensive design of existing Big Data frameworks. It adopts a **task-based** execution engine that is able to **keep CPUs in a cluster fully utilized** even with a budget limit on memory consumption. It also provides a **subgraph-centric API** that is natural for writing graph mining algorithms. G-thinker is orders of magnitude faster than existing systems and scales to graphs orders of magnitude larger given the same hardware environment.

G-thinker works perfectly on all applications except for quasi-clique mining, where we find that some tasks can be so expensive that if we maintain local queue for each thread, head of line blocking will happen leading to low CPU utilization. We thus add global task queue and new task stealing strategy to move "big" tasks forward ASAP.

### Documentation
We maintain detailed documentation in our awesome G-thinker system website: http://www.cs.uab.edu/yanda/gthinker.

It contains detailed steps for deployment and programming which also applies to G-thinker+

However, please replace the code base of G-thinker with the code in this project to compile with G-thinker+. Also note that you may want to implement a new user-defined function Comper::is_bigtask(TaskT * task).

### Contributors
Guimu Guo

Da Yan

### Contact
UAB Data Lab (or YanLab@UAB): https://yanlab19870714.github.io/yanda

Email: yanda@uab.edu
