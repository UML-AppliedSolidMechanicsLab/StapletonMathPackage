# StapletonMathPackage – Maintainer Guide (for Future You)

This doc is the **single source of truth** for how to change the code, test it, and publish a new NuGet version. It assumes you’re using **Visual Studio** on Windows and GitHub Actions is already set up to publish to NuGet when you create a **GitHub Release**.

---

## TL;DR – Release in 10 steps
1. **Pull latest** `main` from GitHub.
2. Make your code changes in `RandomMath/…` (targeting **netstandard2.0**).
3. **Run tests** in VS (Test Explorer) or `dotnet test -c Release`.
4. Update docs if needed (`README.md`).
5. **Bump `<Version>`** in `RandomMath/RandomMath.csproj` (e.g., 1.3.2 → 1.3.3).
6. **Commit & Push** to `main`.
7. **Create tag** in VS (e.g., `v1.3.3`) and **push tags**.
8. On GitHub → **Releases → Draft a new release** → select that tag → **Publish**.
9. Watch **Actions**: “Publish to NuGet” should **Restore → Test → Pack → Push**.
10. Verify on **nuget.org** that the new version appears (indexing may take a few minutes).

> If anything fails, see **Troubleshooting** below.

---

## Local Development

### Prereqs
- Visual Studio 2022+ with .NET 8 SDK installed.
- Repo cloned fresh: `git clone https://github.com/UML-AppliedSolidMechanicsLab/StapletonMathPackage.git`

### Build & Test
- **VS**: Build solution, run tests via **Test Explorer**.
- **CLI** (from repo root):
  ```powershell
  dotnet restore
  dotnet build -c Release
  dotnet test  -c Release
  ```

### Project Targets
- Library targets **`netstandard2.0`** (broad compatibility).
- Tests may target `net8.0`; that’s fine.

---

## Versioning & Packaging

### Where to bump the version
- Edit `RandomMath/RandomMath.csproj` → `<Version>` (Semantic Versioning: `MAJOR.MINOR.PATCH`).

### NuGet metadata (already configured)
- `PackageId` = `StapletonMathPackage`
- `PackageLicenseExpression` = `MIT`
- `PackageReadmeFile` = `README.md`
- README and LICENSE are included in the package.

### Local pack (optional sanity check)
```powershell
dotnet pack .\RandomMath\RandomMath.csproj -c Release -o .\nupkgs
```
You should see `StapletonMathPackage.<version>.nupkg` in `nupkgs/`.

---

## Release & Publish Flow (GitHub Actions)

Publishing is automated when you **Publish a GitHub Release**.

1. **Commit & Push** your changes.
2. **Tag** the commit in VS: `vX.Y.Z` → Push tags.
3. GitHub → **Releases → Draft a new release** → pick that tag → **Auto-generate notes** (or write your own) → **Publish**.
4. Actions runs the workflow: Restore → Test → Pack → Push to NuGet using `NUGET_API_KEY`.

> Alternative: If you prefer publishing **on tag push** (skipping the web release step), change the workflow trigger later to `on: push: tags: ['v*']`.

---

## Updating Dependencies
- **MathNet.Numerics** is referenced in the `.csproj`. To update:
  - In VS: **Manage NuGet Packages** → update **MathNet.Numerics**.
  - Re-run tests; publish a patch/minor release if everything passes.

---

## Do **NOT** Commit Secrets
- Keep secrets (like the NuGet API key) out of the repo.
- The CI uses the GitHub Actions secret **`NUGET_API_KEY`** (Settings → Secrets and variables → Actions).
- If a secret ever slips into git history, follow the **Secret Remediation** steps below.

---

## Troubleshooting

### 1) Pack fails with **README not found**
**Error:** `File not found: RandomMath/README.md`
- The project expects the README one level up. Fix is already in the csproj:
  ```xml
  <PackageReadmeFile>README.md</PackageReadmeFile>
  <ItemGroup>
    <None Include="..\README.md" Pack="true" PackagePath="\" />
    <None Include="..\LICENSE"   Pack="true" PackagePath="\" />
  </ItemGroup>
  ```
- Ensure those files exist at repo root.

### 2) Tests fail only in CI (Linux)
- Ensure the test project is SDK-style and references `Microsoft.NET.Test.Sdk` and your test framework (MSTest/xUnit/NUnit).
- Run the same commands locally:
  ```powershell
  dotnet restore
  dotnet test -c Release --verbosity normal
  ```

### 3) Publish job says **nupkg pushed but version unchanged**
- You forgot to bump `<Version>`.
- Or NuGet is still indexing (refresh the page after a minute).

### 4) Pack/publish complains about wrong TFM
- Confirm library targets **`netstandard2.0`** and you packed the **project** not the solution:
  ```powershell
  dotnet pack .\RandomMath\RandomMath.csproj -c Release -o .\nupkgs
  ```

### 5) “Unauthorized” when pushing to NuGet
- Check repo secret `NUGET_API_KEY` exists and hasn’t expired.
- Regenerate on nuget.org if needed and update the secret.

---

## Secret Remediation (if a key was committed)
1. **Revoke** the key on nuget.org.
2. **Rewrite history** with `git filter-repo` to remove the file/string.
3. **Force-push** the cleaned history.
4. **Re-clone** locally; collaborators must do the same.
5. Add patterns to `.gitignore`:
   ```
   APIKeyNuGet.txt
   *.apikey
   *.secret.txt
   ```

---

## Optional: Accelerated Math Providers
The library ships **managed-only**. If an app wants faster BLAS/LAPACK, it can add a native provider (MKL/OpenBLAS) and call `MathNet.Numerics.Control.UseNative…()` in the **app**, not here.

---

## Housekeeping Checklist (occasionally)
- Run tests and fix warnings.
- Update `README.md` examples if APIs changed.
- Trim unused files and ensure no `Mapack/` remnants are included.
- Consider creating a GitHub **Release** with meaningful notes (helps users and your future self).

---

## Contact & Ownership
- Maintainer: **Scott Stapleton** (UML Applied Solid Mechanics Lab)
- Repo: https://github.com/UML-AppliedSolidMechanicsLab/StapletonMathPackage

Happy hacking! 🔧📦
