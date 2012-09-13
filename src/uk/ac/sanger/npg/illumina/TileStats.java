/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.sanger.npg.illumina;

public class TileStats {
    final int laneNr;
    final int tileNr;
    final int totalCluster;
    final int totalPFCluster;

    public TileStats(int laneNr, int tileNr, int totalCluster, int totalPFCluster) {
        this.laneNr = laneNr;
        this.tileNr = tileNr;
        this.totalCluster = totalCluster;
        this.totalPFCluster = totalPFCluster;
    }
       
    public String toString(){
        return laneNr+"\t"+tileNr+"\t"+totalCluster+"\t"+totalPFCluster;
    }
    
}