
// File: utilities.js
// Version: v1.2
// Date: 2021-02-11
// Authors: Mullissa A., Vollrath A.,  Reiche J., Slagter B., Balling J. , Gou Y., Braun, C.

//---------------------------------------------------------------------------//
// Linear to db scale 
//---------------------------------------------------------------------------//
  
/** Convert backscatter from linear to dB. */
exports.lin_to_db = function(image) {
  var bandNames = image.bandNames().remove('angle');
  var db = ee.Image.constant(10).multiply(image.select(bandNames).log10()).rename(bandNames)
  return image.addBands(db, null, true)
};

/** Convert backscatter from linear to dB. */
exports.db_to_lin = function(image) {
  var bandNames = image.bandNames().remove('angle');
  var lin = ee.Image.constant(10).pow(image.select(bandNames).divide(10)).rename(bandNames)
  return image.addBands(lin, null, true)
};

/*Converts the linear image to db by excluding the ratio bands */
exports.lin_to_db2 = function(image) {
  var db = ee.Image.constant(10).multiply(image.select(['VV', 'VH']).log10()).rename(['VV', 'VH']);
  return image.addBands(db, null, true)
}

//---------------------------------------------------------------------------//
// Prepare ratio band for linear image
//---------------------------------------------------------------------------//
exports.add_ratio_lin = function(image){
      var ratio = image.addBands(image.select('VV').divide(image.select('VH')).rename('VVVH_ratio'));
      return ratio.set('system:time_start', image.get('system:time_start'))
  }


//---------------------------------------------------------------------------//
// Export assets
//---------------------------------------------------------------------------//
var get_options = function(def, options) {
  // default values if undefined
  if (options !== undefined) {
    var opt = options
    for (var key in def) {
      var value = def[key]
      if (opt[key] === undefined) {opt[key] = value}
    }
  } else {var opt = def}
  return opt
}

/*exports individual images to an asset. 
adopted from https://github.com/fitoprincipe/geetools-code-editor */
var Download = {'ImageCollection': {}, 'Table': {}}
Download.ImageCollection.toAsset = function(collection, assetFolder, options) {
  var root = ee.data.getAssetRoots()[0]['id']
  var folder = assetFolder
  if (folder !== null && folder !== undefined) {
    var assetFolder = root+'/'+folder ;
  } else {
    var assetFolder = root
  }
  var defaults = {
      name: null,
      scale: 1000,
      maxPixels: 1e13,
      region: null,
    }
  var opt = get_options(defaults, options)
  var n = collection.size().getInfo();
    
  var colList = collection.toList(n);
  
  var colID = opt.name || collection.getInfo()['id'] || "ImageCollection"
  colID = colID.replace('/','_')
  
  for (var i = 0; i < n; i++) {
    var img = ee.Image(colList.get(i));
    var id = img.id().getInfo() || 'image_'+i.toString();
    var region = opt.region || img.geometry().bounds().getInfo()["coordinates"];
    
    Export.image.toAsset({
      image: img,
      description: id,
      assetId: assetFolder+colID+'_'+id,
      region: region,
      scale: opt.scale,
      maxPixels: opt.maxPixels})
  }
}

exports.Download = Download;

