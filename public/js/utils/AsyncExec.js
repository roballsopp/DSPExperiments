function AsyncExec () {

	var fn = arguments[0];
	var args = [];
	var transferables = [];

	function getTransferables(obj, tfs) {
		for (var prop in obj) {
			if (typeof obj[prop] == 'object') getTransferables(obj[prop], tfs);
			else if (obj[prop] instanceof ArrayBuffer) tfs.push(obj[prop]);
		}
	}

	for (var i = 1, l = arguments.length; i < l; i++) {
		args.push(arguments[i]);
		getTransferables(arguments[i], transferables);
	}

	// Build a worker from an anonymous function body
	var blobURL = URL.createObjectURL( new Blob([ 
		'addEventListener(\'message\', function (e) {\n',
			'console.log("Recieved.");\n',
			'var transferables = [];\n',
			getTransferables.toString() + '\n',
			'var result = ' + fn.toString() + '.apply(this, e.data);\n',
			'getTransferables(result, transferables);\n',
			'console.log("Responding...");\n',
			'postMessage(result, transferables);\n',
		'});'
	], { type: 'application/javascript' }) ),

	worker = new Worker( blobURL );
	worker.postMessage(args, transferables);

	// Won't be needing this anymore
	URL.revokeObjectURL( blobURL );

	return new Promise(function (resolve, reject) {

		worker.addEventListener('message', function (e) {
			resolve(e.data);
			worker.removeEventListener('message', this);
			worker = null;
		});

		worker.addEventListener('error', function (e) {
			reject(e);
			worker.removeEventListener('error', this);
			worker = null;
		});

	});

}
