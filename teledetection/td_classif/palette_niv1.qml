<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.22.4-Białowieża" maxScale="0" styleCategories="AllStyleCategories" minScale="1e+08" hasScaleBasedVisibilityFlag="0">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
    <Private>0</Private>
  </flags>
  <temporal mode="0" fetchMode="0" enabled="0">
    <fixedRange>
      <start></start>
      <end></end>
    </fixedRange>
  </temporal>
  <customproperties>
    <Option type="Map">
      <Option value="false" name="WMSBackgroundLayer" type="bool"/>
      <Option value="false" name="WMSPublishDataSourceUrl" type="bool"/>
      <Option value="0" name="embeddedWidgets/count" type="int"/>
      <Option value="Value" name="identify/format" type="QString"/>
    </Option>
  </customproperties>
  <pipe-data-defined-properties>
    <Option type="Map">
      <Option value="" name="name" type="QString"/>
      <Option name="properties"/>
      <Option value="collection" name="type" type="QString"/>
    </Option>
  </pipe-data-defined-properties>
  <pipe>
    <provider>
      <resampling zoomedOutResamplingMethod="nearestNeighbour" zoomedInResamplingMethod="nearestNeighbour" maxOversampling="2" enabled="false"/>
    </provider>
    <rasterrenderer type="paletted" opacity="1" alphaBand="-1" nodataColor="" band="1">
      <rasterTransparency/>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
        <cumulativeCutLower>0.02</cumulativeCutLower>
        <cumulativeCutUpper>0.98</cumulativeCutUpper>
        <stdDevFactor>2</stdDevFactor>
      </minMaxOrigin>
      <colorPalette>
        <paletteEntry alpha="255" color="#1b7ee1" value="1" label="Eau"/>
        <paletteEntry alpha="255" color="#c91129" value="2" label="Bati"/>
        <paletteEntry alpha="255" color="#ecee44" value="3" label="Sol nu"/>
        <paletteEntry alpha="255" color="#3dc973" value="4" label="Culture"/>
        <paletteEntry alpha="255" color="#0a9534" value="5" label="Ligneux"/>
      </colorPalette>
      <colorramp name="[source]" type="randomcolors">
        <Option/>
      </colorramp>
    </rasterrenderer>
    <brightnesscontrast gamma="1" contrast="0" brightness="0"/>
    <huesaturation colorizeStrength="100" grayscaleMode="0" colorizeBlue="128" colorizeOn="0" invertColors="0" colorizeGreen="128" colorizeRed="255" saturation="0"/>
    <rasterresampler maxOversampling="2"/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
