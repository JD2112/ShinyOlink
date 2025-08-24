# Frequently Asked Questions (FAQ)

??? question "What should I do if the app looks frozen during analysis?"
    A progress bar will always appear while calculations are in progress.  
    If it stops or an error occurs, an error message will be displayed.  

??? question "How should my input files be formatted?"
    You need three files:  
    - **NPX Data (CSV):** Must contain `SampleID`, `Assay`, and `NPX` columns.  
    - **Key File (CSV, optional):** Maps `SampleID` to `SUBJID`.  
    - **Variables File (CSV):** Contains additional variables (e.g., `SUBJID`, Treatment, TimePoint).  

??? question "Can I test the app without my own data?"
    Yes. Download the demo dataset (NPX Data, Key File, Variables File) from  
    [OlinkWrapper on SourceForge](https://sourceforge.net/projects/olinkwrapper/).  

---

## Preprocessing

??? question "Why do I need a Bridge Selector?"
    Bridge proteins ensure consistency across multiple Olink panels.  
    Selecting the right bridge protein is crucial for quality control.  

??? question "What is normalization and why is it important?"
    Normalization adjusts NPX values to reduce technical variation, making samples comparable.  

??? question "How do I handle proteins below the limit of detection (LOD)?"
    Use the **LOD filtering** step to remove proteins with NPX values below the assay’s detection limit.  
    This ensures that only reliable protein measurements are included.  

??? question "How do I detect and handle outliers?"
    The Outlier Detection tool identifies unusual samples or assays.  
    You can remove them if they compromise your analysis.  

---

## Statistical Analysis

??? question "Which test should I use to compare two groups?"
    - Use a **T-test** if your data is normally distributed.  
    - Use a **Wilcoxon test** if your data is non-normal or ordinal.  

??? question "What if I have more than two groups?"
    Use **ANOVA** to test differences across three or more groups.  
    If significant, run a **Post-hoc test** to see which groups differ.  

??? question "What if my samples are not independent (e.g., repeated measures)?"
    Use **Linear Mixed Effects (LME)** models, which handle random effects like subject IDs in longitudinal studies.  

---

## Visualization

??? question "Why is my Volcano Plot empty?"
    You must first run a **T-test** or **Wilcoxon test**.  
    The Volcano Plot depends on those results.  

??? question "Why can’t I generate a Pathway Heatmap?"
    You must run a **Pathway Enrichment Analysis** first.  
    The heatmap relies on enrichment results.  

??? question "How can I explore global protein patterns?"
    Use **PCA** or **UMAP** plots to visualize high-dimensional data in 2D/3D.  

??? question "What plots are available for presenting results?"
    The app supports:  
    - Box plots  
    - Violin plots  
    - Heatmaps  
    - QC plots  
    - LME plots  
    - Volcano plots  
    - Pathway heatmaps  
    - Distribution plots  

---

## Advanced Analysis

??? question "What is Pathway Enrichment Analysis?"
    It tests whether your significant proteins (from T-test or Wilcoxon) are overrepresented in known biological pathways.  

??? question "Can I run linear regression in the app?"
    Yes. Linear regression models protein NPX levels against continuous variables (e.g., age, BMI, clinical scores).  

---

## Server & Browser Issues

??? question "I see a 404 error when accessing the app. What should I do?"
    - Wait a few moments and refresh the browser.  
    - If the problem persists, clear your browser cache and cookies.  
    - Restart your browser or computer before retrying.  
    - If you still cannot access the app, contact your system administrator.  

??? question "The app is very slow or not responding. How can I fix this?"
    - Ensure you have a stable internet connection.  
    - Try using a different browser (Google Chrome or Firefox recommended).  
    - Close unnecessary browser tabs or applications consuming high memory.  
    - If running locally, check that your machine has sufficient RAM/CPU for the dataset size.  

??? question "The plots or tables are not displaying properly."
    - Refresh the browser page.  
    - Clear the browser cache and reload.  
    - Try resizing the browser window or switching to full-screen mode.  
    - If the problem continues, restart the app/server.  

??? question "The app logged me out unexpectedly."
    - This can happen if your session times out.  
    - Refresh the browser and log in again.  
    - If the issue repeats frequently, ask your administrator to extend the session timeout.
    - It can also occur due to the ill-formatted input files. Please double-check your files with the required format (See "Demo Data" section).

??? question "I uploaded my files but nothing happens."
    - Double-check that the files are in the correct **CSV format**.  
    - Ensure your **NPX Data**, **Key File**, and **Variables File** follow the required structure.  
    - Re-upload the files one by one and then click **Merge Data**.  
    - If the issue persists, try using a different browser or clear your current browser's cache.  
    - Contact support if you continue to experience issues.