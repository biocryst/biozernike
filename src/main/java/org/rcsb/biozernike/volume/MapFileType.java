package org.rcsb.biozernike.volume;

public enum MapFileType {

    CCP4(".map"),
    MRC(".mrc");

    private final String extension;

    private MapFileType(String extension) {
        this.extension = extension;
    }

    public String getExtension() {
        return extension;
    }
}
