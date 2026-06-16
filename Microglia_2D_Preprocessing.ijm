// ============================================================
// Microglia 2D Preprocessing Macro
// Automates the ImageJ portion of the Microglia 2D Sholl
// Analysis SOP (A. Schwarz, Version 1.0, 2/19/26)
//
// Steps performed per .czi file:
//   1. Open with Bio-Formats (split channels, autoscale)
//   2. Discard PNN and PV channels
//   3. Subtract background (rolling ball) from glia and DAPI
//   4. Z-project both channels (max intensity)
//   5. Merge into a 2-channel composite (C1=glia, C3=DAPI)
//   6. Save as TIFF
// ============================================================

// ---- Configuration ----
// Channel indices (0-based) in the source .czi files
// Modify these if your acquisition channel order differs
glia_channel = 0;    // IBA1 / microglia
dapi_channel = 2;    // DAPI nuclear stain
// Channels at indices 1 (PNN) and 3 (PV) will be discarded

rolling_ball_radius = 100;  // pixels, for background subtraction

// ---- Directory Selection ----
input_dir = getDirectory("Select folder containing .czi files");
output_dir = getDirectory("Select folder for saving 2D composites");

// ---- Batch Mode ----
setBatchMode(true);

// ---- Get File List ----
file_list = getFileList(input_dir);

// ---- Process Each File ----
for (i = 0; i < file_list.length; i++) {

    // Skip any file that isn't a .czi
    if (!endsWith(file_list[i], ".czi")) {
        continue;
    }

    filename = file_list[i];
    basename = replace(filename, ".czi", "");

    print("Processing: " + filename);

    // Open .czi with Bio-Formats, splitting channels
    run("Bio-Formats Importer",
        "open=[" + input_dir + filename + "] " +
        "autoscale " +
        "color_mode=Default " +
        "split_channels " +
        "view=Hyperstack " +
        "stack_order=XYCZT");

    // ---- Discover actual window titles ----
    // Bio-Formats may name windows differently depending on
    // version and file structure, so we find them by matching
    // the " - C=N" suffix rather than constructing names.
    titles = getList("image.titles");
    glia_title = "";
    dapi_title = "";

    for (t = 0; t < titles.length; t++) {
        if (endsWith(titles[t], " - C=" + glia_channel)) {
            glia_title = titles[t];
        }
        if (endsWith(titles[t], " - C=" + dapi_channel)) {
            dapi_title = titles[t];
        }
    }

    // Verify we found both channels
    if (glia_title == "" || dapi_title == "") {
        print("ERROR: Could not find expected channels in " + filename);
        print("  Looking for channels C=" + glia_channel + " and C=" + dapi_channel);
        print("  Windows found after import:");
        for (t = 0; t < titles.length; t++) {
            print("    " + titles[t]);
        }
        // Close all windows from this file before skipping
        while (nImages > 0) { close(); }
        continue;
    }

    // Close all channels we are NOT keeping
    for (t = 0; t < titles.length; t++) {
        if (titles[t] != glia_title && titles[t] != dapi_title) {
            close(titles[t]);
        }
    }

    // ---- Background Subtraction ----
    // Process glia (IBA1) channel
    selectWindow(glia_title);
    run("Subtract Background...", "rolling=" + rolling_ball_radius + " stack");

    // Process DAPI channel
    selectWindow(dapi_title);
    run("Subtract Background...", "rolling=" + rolling_ball_radius + " stack");

    // ---- Z Projection ----
    // Project glia channel
    selectWindow(glia_title);
    n_slices = nSlices;
    run("Z Project...", "start=1 stop=" + n_slices + " projection=[Max Intensity]");
    glia_projected = getTitle();
    close(glia_title);

    // Project DAPI channel
    selectWindow(dapi_title);
    n_slices = nSlices;
    run("Z Project...", "start=1 stop=" + n_slices + " projection=[Max Intensity]");
    dapi_projected = getTitle();
    close(dapi_title);

    // ---- Merge Channels ----
    run("Merge Channels...",
        "c1=[" + glia_projected + "] " +
        "c3=[" + dapi_projected + "] " +
        "create");
    composite_title = getTitle();

    // ---- Save Result ----
    output_filename = basename + "-2D-IBA1-DAPI.tiff";
    saveAs("Tiff", output_dir + output_filename);
    print("Saved: " + output_filename);

    // Close the saved image before processing the next file
    close();
}

// ---- Cleanup ----
setBatchMode(false);
print("---- Batch complete ----");
print("Processed files saved to: " + output_dir);
