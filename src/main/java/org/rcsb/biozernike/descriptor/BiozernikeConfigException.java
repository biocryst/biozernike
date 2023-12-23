package org.rcsb.biozernike.descriptor;

public class BiozernikeConfigException extends RuntimeException {

    private static final long serialVersionUID = -4278211404060006193L;

    public BiozernikeConfigException(String message) {
        super(message);
    }

    public BiozernikeConfigException(String message, Throwable cause) {
        super(message, cause);
    }
}
