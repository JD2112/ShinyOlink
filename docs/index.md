# Tutorial: ShinyOlink

ShinyOlink is a user-friendly **Shiny application** designed to simplify the *analysis of Olink proteomics data*. It allows users to **upload their raw data**, then performs various *statistical tests* and generates *interactive visualizations* by leveraging the specialized `OlinkAnalyze` R package, making complex bioinformatics accessible to researchers. The project uses a *modular design* and has *automated deployment pipelines* for efficient delivery and maintenance.


## Visual Overview

```mermaid
flowchart TD
    A0["Shiny Application Structure"]
    A1["Module-Based Development"]
    A2["Data Ingestion & Preparation"]
    A3["OlinkAnalyze Package Integration"]
    A4["Reactive Programming"]
    A5["Statistical Analysis & Visualization Modules"]
    A6["Automated Deployment Pipelines"]
    A7["Project Governance & Contributions"]
    A0 -- "Composed of Modules" --> A1
    A0 -- "Underpins Application Logic" --> A4
    A1 -- "Implements Functionality" --> A5
    A1 -- "Manages Data Upload" --> A2
    A2 -- "Provides Reactive Data" --> A4
    A2 -- "Supplies Data For" --> A5
    A3 -- "Provides Analysis Functions" --> A5
    A4 -- "Enables Dynamic Outputs" --> A5
    A6 -- "Automates Deployment" --> A0
    A7 -- "Governs Deployment" --> A6
    A7 -- "Guides Module Design" --> A1
```

## Chapters

1. [Data Ingestion & Preparation](chapter1.md)
2. [Statistical Analysis & Visualization Modules](chapter2.md)
3. [OlinkAnalyze Package Integration](chapter3.md)
4. [Shiny Application Structure](chapter4.md)
5. [Module-Based Development](chapter5.md)
6. [Reactive Programming](chapter6.md)
7. [Automated Deployment Pipelines](chapter7.md)
8. [Project Governance & Contributions](chapter8.md)

---