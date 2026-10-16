# The following code is the most recent version of combined Hydrocode <u>RocfluMP</u> and particle-tracker <u>ppiclF</u>.

## 📘 **Documentations**: 
### RocfluMP
You can find the full RocfluMP manual [here](docs/Rocflu_manual.pdf).

### ppiclF
Documentation Website for ppiclF [here](https://dpzwick.github.io/ppiclF-doc/user/external.html).

## Getting Started:

## 🚀 Getting Started

Follow these steps to get a copy of the project up and running on your machine or HiPerGator account.

### 1. 🍴 Fork the Repository

- Click the **Fork** button at the top-right of this page to create your own copy of the repository.

### 2. 📥 Clone Your Fork

#### Option A: Clone via **HTTPS**
```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```

#### Option B: Clone via **SSH** (recommended for contributors)

1. **Generate an SSH key** (if you don’t already have one):

   ```bash
   ssh-keygen -t ed25519 -C "your_email@example.com"
   ```

2. **Add your SSH key to the ssh-agent**:

   ```bash
   eval "$(ssh-agent -s)"
   ssh-add ~/.ssh/id_ed25519
   ```

3. **Add the SSH key to your GitHub account**:
   Display your public key with:
   ```bash
   cat ~/.ssh/id_ed25519.pub
   ```
   Go to GitHub → Settings → SSH and GPG keys → New SSH key
   
   Paste the copied key into the field and give it a descriptive title (e.g., My Laptop or HiPerGator).

4. **Clone your fork using SSH**:
   ```bash
   git clone git@github.com:your-username/your-repo-name.git
   cd your-repo-name
   ```

## 🔨 Building

The code builds with CMake (3.20 or newer). On HiPerGator:

```bash
module purge; module load gcc/14.2.0 openmpi/5.0.7
cmake -S . -B build
cmake --build build -j 16
```

The executables go to `build/bin`: the solver `rflump` and the three
pre-processing utilities `rflumap`, `rflupart` and `rfluinit` that set up a case.
A clean build takes about 3 minutes with 16 cores. After an edit, run the
`cmake --build` line again: only the changed files and the files that depend on
them are recompiled.

Options are set with `-D` on the first `cmake` line:

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Debug`: `-g -O0 -fbacktrace -fcheck=all` everywhere |
| `ROCFLU_PICL` | `ON` | couple the ppiclF particle solver |
| `ROCFLU_SPEC` | `OFF` | species transport (also needed by the JWL gas model, `GASMODEL 7`) |
| `ROCFLU_PLAG` | `OFF` | Rocpart Lagrangian particles |
| `ROCFLU_FOLDER` | `ON` | write output into per-kind folders |
| `ROCFLU_UTILITIES` | `ON` | build `rflumap`, `rflupart` and `rfluinit` |
| `METIS_ROOT` | `3rd_party/local` | METIS 5 and GKlib install (`3rd_party/install.sh` builds them) |

For example, a debug build with species in its own directory:
`cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DROCFLU_SPEC=ON`.

Release builds compile Rocflu with `-O` and ppiclF with `-O3 -fcheck=all`.
The sources each part compiles are listed in the `Filelist.txt` of its
directory. The old make build (`makeall.sh`) still works; its executables give
results identical to the CMake ones, byte for byte.

## 🧪 Testing

The regression suite in [`tests/`](tests/README.md) runs 35 small cases of the coupled
code on 1, 2 and 4 MPI ranks. It compares each one with a stored result and checks the
physics (conservation, interpolation, collisions, walls). After building:

```bash
cd tests
module purge; module load gcc/14.2.0 openmpi/5.0.7
./run_tests.py run --exe ../build/bin
```

See [`tests/README.md`](tests/README.md) for the list of cases and how to run them on
HiPerGator, and [`docs/PseudoTurb_fixes.md`](docs/PseudoTurb_fixes.md) for the fixes
they guard.

## 🛠️ Making Changes with Git

Once you’ve cloned your fork of the repository, here are the basic Git commands you’ll use when contributing to the project:

### 🔄 1. Pull the Latest Changes

Before making changes, make sure your local branch is up to date with the main project:

```bash
git fetch upstream
git checkout main
git merge upstream/main
```

### 📝 2. Create a New Branch for Your Changes
It's good practice to make changes on a separate branch, not directly on `main`.
```bash
git checkout -b feature/your-branch-name
```

### 💻 3. Make Your Code Changes
Edit the files you need using your favorite editor or IDE.

### ✅ 4. Stage and Commit Your Changes
```bash 
git add .
git commit -m "Descriptive message of what you changed"
```

### 📤 5. Push to Your Fork
```branch
git push origin feature/your-branch-name
```

### 🚀 6. Open a Pull Request
1. Go to your forked repository on GitHub.

2. You’ll see a prompt to open a pull request (PR) from your branch.

3. Add a clear description of your changes and submit.


## 🧼 Recommended Tips


- Use `.gitignore` (already found in your cloned directory) to avoid committing unwanted files.

- Run `git status` often to check what’s staged or modified.

- Run `git log --oneline` to view your commit history.
