package es.uam.eps.ir.sr.utils;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


/***
 * 
 */

public class FrequencyGrouping {
    private List<Double> lowGroup;
    private List<Double> midGroup;
    private List<Double> highGroup;
    private double lowGroupPercentage;
    private double midGroupPercentage;
    private double highGroupPercentage;
    
    private double lowGroupMedian;
    private double midGroupMedian;
    private double highGroupMedian;
    
    private int maxSize;

    public FrequencyGrouping(int maxSize) {
        reset();
        this.maxSize = maxSize;
    }
    
    private void reset() {
    	this.lowGroup = new ArrayList<>();
        this.midGroup = new ArrayList<>();
        this.highGroup = new ArrayList<>();
        this.lowGroupPercentage = 0.0;
        this.midGroupPercentage = 0.0;
        this.highGroupPercentage = 0.0;
        
        this.lowGroupMedian = 0;
        this.midGroupMedian = 0;
        this.highGroupMedian = 0;
    }
    
    public double giveValue(int counter) {
    	
    	double actualpercent = (counter * 100) / (this.maxSize);
    	
    	if (actualpercent <= this.highGroupPercentage) {
    		return this.highGroupMedian;
    	}
    	else if (actualpercent <= this.highGroupPercentage + this.midGroupPercentage) {
    		return this.midGroupMedian;
    	}
    	else
    		return this.lowGroupMedian;
    }
    

    public void groupFrequenciesBySize(List<Double> frequencies) {
    	reset();
    	
        if (frequencies == null || frequencies.isEmpty()) {
            throw new IllegalArgumentException("La lista de frecuencias no puede estar vac�a.");
        }

        // Sort
        Collections.sort(frequencies);

        // Thresholds
        double lowLimit = frequencies.get(frequencies.size() / 3);
        double highLimit = frequencies.get(2 * frequencies.size() / 3);

        // Grouping
        for (double frequency : frequencies) {
            if (frequency <= lowLimit) {
                lowGroup.add(frequency);
            } else if (frequency <= highLimit) {
                midGroup.add(frequency);
            } else {
                highGroup.add(frequency);
            }
        }

        int total = frequencies.size();
        lowGroupPercentage = (double) lowGroup.size() / total * 100;
        midGroupPercentage = (double) midGroup.size() / total * 100;
        highGroupPercentage = (double) highGroup.size() / total * 100;
        
        lowGroupMedian = computeMedian(lowGroup);
        midGroupMedian = computeMedian(midGroup);
        highGroupMedian = computeMedian(highGroup);
    }
    
    /**
    private int computeMean(List<Double> group) {
        if (group.isEmpty()) {
            return 0;
        }
        double sum = 0.0;
        for (double value : group) {
            sum += value;
        }
        return (int) Math.round(sum / group.size());    
    }
    **/
    private double computeMedian(List<Double> group) {
    	if (group.isEmpty()) {
            return 1.0; //I use 1.0 because all the files contains checckins, so the default value should be 1.0
        }
    	
        Collections.sort(group);
        
        int tam = group.size();
        double median;

        if (tam % 2 == 0) {
        	median = (group.get(tam / 2 - 1) + group.get(tam / 2)) / 2.0;
        } else {
        	median = group.get(tam / 2);
        }

        return median;
    }
 

    
    
    public void groupFrequenciesByValues(List<Double> frequencies) {
    	reset();
        
        if (frequencies == null || frequencies.isEmpty()) {
            throw new IllegalArgumentException("La lista de frecuencias no puede estar vac�a.");
        }

        // Sort
        Collections.sort(frequencies);

        // Thresholds
        double highThreshold = frequencies.get(frequencies.size() - 1);
        double lowThreshold = frequencies.get(0);

        // Limit of grouping
        double range = highThreshold - lowThreshold;
        double midCutoff = lowThreshold + range / 3;
        double highCutoff = lowThreshold + 2 * range / 3;

        // Agrupar las frecuencias en "bajos", "medios" y "altos"
        for (double frequency : frequencies) {
            if (frequency > highCutoff) {
                highGroup.add(frequency);
            } else if (frequency > midCutoff) {
                midGroup.add(frequency);
            } else {
                lowGroup.add(frequency);
            }
        }

        // Calcular los porcentajes
        int total = frequencies.size();
        lowGroupPercentage = (double) lowGroup.size() / total * 100;
        midGroupPercentage = (double) midGroup.size() / total * 100;
        highGroupPercentage = (double) highGroup.size() / total * 100;
        
        lowGroupMedian = computeMedian(lowGroup);
        midGroupMedian = computeMedian(midGroup);
        highGroupMedian = computeMedian(highGroup);
    }

    public List<Double> getLowGroup() {
        return lowGroup;
    }

    public List<Double> getMidGroup() {
        return midGroup;
    }

    public List<Double> getHighGroup() {
        return highGroup;
    }

    public int getLowGroupSize() {
        return lowGroup.size();
    }

    public int getMidGroupSize() {
        return midGroup.size();
    }

    public int getHighGroupSize() {
        return highGroup.size();
    }

    public double getLowGroupPercentage() {
        return lowGroupPercentage;
    }

    public double getMidGroupPercentage() {
        return midGroupPercentage;
    }

    public double getHighGroupPercentage() {
        return highGroupPercentage;
    }
    
    public double getLowGroupMedian() {
        return lowGroupMedian;
    }

    public double getMidGroupMedian() {
        return midGroupMedian;
    }

    public double getHighGroupMedian() {
        return highGroupMedian;
    }

    public static void main(String[] args) {
        FrequencyGrouping frequencyGrouping = new FrequencyGrouping(5);
        
        
        List<Double> frequencies = new ArrayList<>(List.of(10.5, 23.4, 56.7, 12.3, 67.8, 45.6, 34.5, 78.9, 21.2, 43.1));
        
        frequencyGrouping.groupFrequenciesBySize(frequencies);

        System.out.println("Grupo Bajo: " + frequencyGrouping.getLowGroup());
        System.out.println("Grupo Medio: " + frequencyGrouping.getMidGroup());
        System.out.println("Grupo Alto: " + frequencyGrouping.getHighGroup());

        System.out.println("Tama�o del Grupo Bajo: " + frequencyGrouping.getLowGroupSize());
        System.out.println("Tama�o del Grupo Medio: " + frequencyGrouping.getMidGroupSize());
        System.out.println("Tama�o del Grupo Alto: " + frequencyGrouping.getHighGroupSize());

        System.out.println("Porcentaje del Grupo Bajo: " + frequencyGrouping.getLowGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Medio: " + frequencyGrouping.getMidGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Alto: " + frequencyGrouping.getHighGroupPercentage() + "%");
        
        for (int i = 1; i <= 5; i++) {
        	System.out.println(frequencyGrouping.giveValue(i));
        }
        
        
        
        
        frequencyGrouping.groupFrequenciesByValues(frequencies);

        System.out.println("Grupo Bajo: " + frequencyGrouping.getLowGroup());
        System.out.println("Grupo Medio: " + frequencyGrouping.getMidGroup());
        System.out.println("Grupo Alto: " + frequencyGrouping.getHighGroup());

        System.out.println("Tama�o del Grupo Bajo: " + frequencyGrouping.getLowGroupSize());
        System.out.println("Tama�o del Grupo Medio: " + frequencyGrouping.getMidGroupSize());
        System.out.println("Tama�o del Grupo Alto: " + frequencyGrouping.getHighGroupSize());

        System.out.println("Porcentaje del Grupo Bajo: " + frequencyGrouping.getLowGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Medio: " + frequencyGrouping.getMidGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Alto: " + frequencyGrouping.getHighGroupPercentage() + "%");
        
        System.out.println("Media del Grupo Alto: " + frequencyGrouping.getHighGroupMedian());
        System.out.println("Media del Grupo Medio: " + frequencyGrouping.getMidGroupMedian());
        System.out.println("Media del Grupo Bajo: " + frequencyGrouping.getLowGroupMedian());
        
        frequencies = new ArrayList<>(List.of(100.0, 90.0, 80.0, 79.0, 78.0, 20.0, 1.0));

        frequencyGrouping.groupFrequenciesBySize(frequencies);

        System.out.println("Grupo Alto: " + frequencyGrouping.getHighGroup());
        System.out.println("Grupo Medio: " + frequencyGrouping.getMidGroup());
        System.out.println("Grupo Bajo: " + frequencyGrouping.getLowGroup());

        System.out.println("Tama�o del Grupo Alto: " + frequencyGrouping.getHighGroupSize());
        System.out.println("Tama�o del Grupo Medio: " + frequencyGrouping.getMidGroupSize());
        System.out.println("Tama�o del Grupo Bajo: " + frequencyGrouping.getLowGroupSize());

        System.out.println("Porcentaje del Grupo Alto: " + frequencyGrouping.getHighGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Medio: " + frequencyGrouping.getMidGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Bajo: " + frequencyGrouping.getLowGroupPercentage() + "%");
        
        frequencyGrouping.groupFrequenciesByValues(frequencies);

        System.out.println("Grupo Alto: " + frequencyGrouping.getHighGroup());
        System.out.println("Grupo Medio: " + frequencyGrouping.getMidGroup());
        System.out.println("Grupo Bajo: " + frequencyGrouping.getLowGroup());

        System.out.println("Tama�o del Grupo Alto: " + frequencyGrouping.getHighGroupSize());
        System.out.println("Tama�o del Grupo Medio: " + frequencyGrouping.getMidGroupSize());
        System.out.println("Tama�o del Grupo Bajo: " + frequencyGrouping.getLowGroupSize());

        System.out.println("Porcentaje del Grupo Alto: " + frequencyGrouping.getHighGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Medio: " + frequencyGrouping.getMidGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Bajo: " + frequencyGrouping.getLowGroupPercentage() + "%");
        
        System.out.println("Media del Grupo Alto: " + frequencyGrouping.getHighGroupMedian());
        System.out.println("Media del Grupo Medio: " + frequencyGrouping.getMidGroupMedian());
        System.out.println("Media del Grupo Bajo: " + frequencyGrouping.getLowGroupMedian());
        
        frequencies = new ArrayList<>(List.of(1.0));
        
        frequencyGrouping.groupFrequenciesBySize(frequencies);

        System.out.println("Grupo Alto: " + frequencyGrouping.getHighGroup());
        System.out.println("Grupo Medio: " + frequencyGrouping.getMidGroup());
        System.out.println("Grupo Bajo: " + frequencyGrouping.getLowGroup());

        System.out.println("Tama�o del Grupo Alto: " + frequencyGrouping.getHighGroupSize());
        System.out.println("Tama�o del Grupo Medio: " + frequencyGrouping.getMidGroupSize());
        System.out.println("Tama�o del Grupo Bajo: " + frequencyGrouping.getLowGroupSize());

        System.out.println("Porcentaje del Grupo Alto: " + frequencyGrouping.getHighGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Medio: " + frequencyGrouping.getMidGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Bajo: " + frequencyGrouping.getLowGroupPercentage() + "%");
        
        System.out.println("Media del Grupo Alto: " + frequencyGrouping.getHighGroupMedian());
        System.out.println("Media del Grupo Medio: " + frequencyGrouping.getMidGroupMedian());
        System.out.println("Media del Grupo Bajo: " + frequencyGrouping.getLowGroupMedian());
        
        frequencyGrouping.groupFrequenciesByValues(frequencies);

        System.out.println("Grupo Alto: " + frequencyGrouping.getHighGroup());
        System.out.println("Grupo Medio: " + frequencyGrouping.getMidGroup());
        System.out.println("Grupo Bajo: " + frequencyGrouping.getLowGroup());

        System.out.println("Tama�o del Grupo Alto: " + frequencyGrouping.getHighGroupSize());
        System.out.println("Tama�o del Grupo Medio: " + frequencyGrouping.getMidGroupSize());
        System.out.println("Tama�o del Grupo Bajo: " + frequencyGrouping.getLowGroupSize());

        System.out.println("Porcentaje del Grupo Alto: " + frequencyGrouping.getHighGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Medio: " + frequencyGrouping.getMidGroupPercentage() + "%");
        System.out.println("Porcentaje del Grupo Bajo: " + frequencyGrouping.getLowGroupPercentage() + "%");
        
        System.out.println("Media del Grupo Alto: " + frequencyGrouping.getHighGroupMedian());
        System.out.println("Media del Grupo Medio: " + frequencyGrouping.getMidGroupMedian());
        System.out.println("Media del Grupo Bajo: " + frequencyGrouping.getLowGroupMedian());

    }
}